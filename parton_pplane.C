//-----------------------------------------------
//Code for Parton participant plane calculation
//AMPT model
//
//
//Author: P. Yin
//-----------------------------------------------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <iterator>
#include <cmath>
#include <TLegend.h>
#include <TLatex.h>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TStyle.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"

using namespace std;

//-----------------------------------------------------------------------------------
//Structure declarations
struct particle
{
	int   id;
	float px;
	float py;
	float pz;
	float m;
	float x;
	float y;
	float z;
	float t;
	float eta;
	float phi;
	float pT;
	float rsquare;
};

//-----------------------------------------------------------------------------------
//Global variables declarations
int event_counter = 0;

//-----------------------------------------------------------------------------------
//Vector declarations
vector<particle> partons;

vector<float> psi2;

vector<particle> final_p;

//-----------------------------------------------------------------------------------
//Graphs declarations
TProfile* v2s;
TProfile* v2s_pt;

//-----------------------------------------------------------------------------------
//Output file

//-----------------------------------------------------------------------------------
//Test cout turn on/off set up
bool test = false;

//-----------------------------------------------------------------------------------
//Functions declarations
//Calculate parton participant plane system
void  processEvent_melting()
{
	float cmx = 0;
	float cmy = 0;
	float qx2 = 0;
	float qx3 = 0;
	float qy2 = 0;
	float qy3 = 0;
	float aver2 = 0;

	for (unsigned int i = 0; i < partons.size(); i++)
	{
		cmx = cmx + partons[i].x;
		cmy = cmy + partons[i].y;
	}

	int count = partons.size();

	cmx = cmx / (float)count;
	cmy = cmy / (float)count;

	for (unsigned int i = 0; i < partons.size(); i++)
	{
		//Shift to center of mass frame
		partons[i].x = partons[i].x - cmx;
		partons[i].y = partons[i].y - cmy;

		partons[i].phi = TMath::ATan2(partons[i].y, partons[i].x);
		partons[i].rsquare = partons[i].x * partons[i].x + partons[i].y * partons[i].y;
	}

	for (unsigned int i = 0; i < partons.size(); i++)
	{
		qx2 = qx2 + partons[i].rsquare * TMath::Cos(2 * partons[i].phi);
		qy2 = qy2 + partons[i].rsquare * TMath::Sin(2 * partons[i].phi);

		qx3 = qx3 + partons[i].rsquare * TMath::Cos(3 * partons[i].phi);
		qy3 = qy3 + partons[i].rsquare * TMath::Sin(3 * partons[i].phi);

		aver2    = aver2 + partons[i].rsquare;
	}

	qx2 = qx2 / (float)count;
	qy2 = qy2 / (float)count;

	qx3 = qx3 / (float)count;
	qy3 = qy3 / (float)count;

	aver2    = aver2 / (float)count;

	float temp_psi2 = (TMath::ATan2(qy2, qx2) + TMath::Pi()) / 2;

	psi2.push_back(temp_psi2);
}

//parse initial state particles file
void parse_afterPropagation_file()
{
	//Read in data file
	ifstream partonFile;
	//partonFile.open("parton-initial-afterPropagation.dat");
	partonFile.open("ana/parton-initial-afterPropagation.dat");

	//Skip the job if not partonFile
	if (!partonFile)
	{
		cout << "-->No such data file" << endl << endl;
		return;
	}
	else cout << "--> Successfully open the afterPropagation file!" << endl << endl;

	//In this while loop, program will read the data file line by line
	while (partonFile)
	{
		int    evtnumber;
		int    iteration;
		int    nlist;
		int    n_baryon_formed;
		int    n_meson_formed;
		int    n_inipart;
		int    n_inipart_notinzpc;

		//Get the header of each event
		partonFile >> evtnumber >> iteration >> nlist >> n_baryon_formed >> n_meson_formed >> n_inipart >> n_inipart_notinzpc;

		if (!partonFile) break;

		//Analysis each particle in the event
		for (int i = 0; i < nlist; i++)
		{
			int    partid;
			float  pv[3];
			float  mass;
			double space[4];

			partonFile >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

			if (TMath::Sqrt(pv[0] * pv[0] + pv[1] * pv[1]) < 0.0001) continue;

			//Calculate the energy
			float energy = TMath::Sqrt(pv[0] * pv[0] + pv[1] * pv[1] + pv[2] * pv[2] + mass * mass);

			//Make Lorentz vector
			TLorentzVector ev(pv[0], pv[1], pv[2], energy);

			//Get pT, phi, pseudorapidity, particle id, px, py, and pz. Store them into p.
			particle p;

			p.eta = ev.Eta();
			p.pT  = ev.Pt();
			p.px  = pv[0];
			p.py  = pv[1];
			p.pz  = pv[2];
			p.m   = mass;
			p.x   = space[0];
			p.y   = space[1];
			p.z   = space[2];
			p.t   = space[3];

			if (abs(p.eta) < 3 && p.t < 3) partons.push_back(p);
		}

		processEvent_melting();

		partons.clear();
	}
}

//Calculate v2 vs. pt and N_charge
void  processEvent_ampt(int evtnumber, int ncharge)
{
	//Calculate v2
	for (unsigned int i=0; i<final_p.size(); i++)
	{
		float v2 = TMath::Cos(2 * (final_p[i].phi - psi2[evtnumber]));
		v2s->Fill(ncharge, v2);
		v2s_pt->Fill(final_p[i].pT, v2);
	}
}

//parse final state partilces code
void parse_ampt_file()
{
    //Read in data file
    ifstream dataFile;
    //dataFile.open("ampt.dat");
    dataFile.open("ana/ampt.dat");

    //Skip the job if not dataFile
    if (!dataFile)
    {
        cout << "-->No such data file" << endl << endl;
        return;
    }
    else cout << "--> Successfully open the data file!" << endl << endl;

    //In this while loop, program will read the data file line by line
    while (dataFile)
    {
        int    evtnumber;
        int    testnum;
        int    nlist;
        double impactpar;
        int    npartproj;
        int    nparttarg;
        int    npartprojelas;
        int    npartprojinelas;
        int    nparttargelas;
        int    nparttarginelas;
        double junk;

        int    ncharge = 0;

        //Get the header of each event
        dataFile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk;

        if (!dataFile) break;

        //Analysis each particle in the event
        for (int i = 0; i < nlist; i++)
        {
            int    partid;
            float  pv[3];
            float  mass;
            double space[4];

            dataFile >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

            //Skip non-charged particles that we are not interested
            //
            //+-211 --> pions,  +-321 --> kaons, +-2212 --> protons
            if (abs(partid) != 211 && abs(partid) != 321 && abs(partid) != 2212) continue;

            if (TMath::Sqrt(pv[0] * pv[0] + pv[1] * pv[1]) < 0.0001) continue;

            //Calculate the energy
            float energy = TMath::Sqrt(pv[0] * pv[0] + pv[1] * pv[1] + pv[2] * pv[2] + mass * mass);

            //Make Lorentz vector
            TLorentzVector ev(pv[0], pv[1], pv[2], energy);

            //Get pT, phi, pseudorapidity, particle id, px, py, and pz. Store them into p.
            particle p;

            p.eta = ev.Eta();
            p.pT  = ev.Pt();
            p.phi = ev.Phi();
            p.px  = pv[0];
            p.py  = pv[1];
            p.pz  = pv[2];

            if (abs(p.eta) > 1 && abs(p.eta) < 3  && p.pT > 0.3 && p.pT < 3) 
            {
            	final_p.push_back(p);
            	ncharge++;
            }

        }
        processEvent_ampt(event_counter, ncharge);

        final_p.clear();

        event_counter += 1;

        if (!dataFile) break;
    }
}


void parton_pplane()
{
    TFile *fout = new TFile("parton_pplane.root", "RECREATE");

    // 60, -0.5, 599.5, -10, 10
    // 50, -0.5, 499.5, -10, 10     
    // 100, -0.5, 5999.5, -10, 10 Pb+Pb
    // 50, -0.5, 199.5, -10, 10    d+Au
    v2s = new TProfile("v2s", "v2s", 50, -0.5, 199.5, -10, 10); //
    v2s_pt = new TProfile("v2s_pt", "v2s_pt", 13, 0.2, 3, -1.0, 1.0);

    parse_afterPropagation_file();
    parse_ampt_file();

    fout->Write();
    fout->Close();
}













































