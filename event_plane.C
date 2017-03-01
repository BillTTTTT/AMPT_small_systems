//-----------------------------------------------
//Code to run AMPT
//event plane method
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

//-----------------------------------------------------------------------------------
//Vector declarations
vector<particle> total_particles;
//BBC [-3.9,-3)
vector<particle> pBBCS;

//FVTXS [-3,-1)
vector<particle> pFVTXS;

//FVTXN [1,3)
vector<particle> pFVTXN;

//-----------------------------------------------------------------------------------
//Graphs declarations
TProfile* v2s; //v2s
TProfile* res_comp;

//-----------------------------------------------------------------------------------
//Output file

//-----------------------------------------------------------------------------------
//Test cout turn on/off set up
bool test = false;

//-----------------------------------------------------------------------------------
//Functions declarations
void  processEvent(vector<particle> pA, vector<particle> pB, vector<particle> pC)
{
    if (pA.size() == 0) return;

    float qxA = 0;
    float qyA = 0;
    float qxB = 0;
    float qyB = 0;
    float qxC = 0;
    float qyC = 0;
    float psiA = 0;
    float psiB = 0;
    float psiC = 0;

    for (unsigned int i = 0; i < pA.size(); i++)
    {
        qxA += TMath::Cos(2 * pA[i].phi);
        qyA += TMath::Sin(2 * pA[i].phi);
    }

    qxA = qxA / pA.size();
    qyA = qyA / pA.size();
    psiA = TMath::ATan2(qyA, qxA) / 2;

    for (unsigned int i = 0; i < pB.size(); i++)
    {
        qxB += TMath::Cos(2 * pB[i].phi);
        qyB += TMath::Sin(2 * pB[i].phi);
    }

    qxB = qxB / pB.size();
    qyB = qyB / pB.size();
    psiB = TMath::ATan2(qyB, qxB) / 2;

    for (unsigned int i = 0; i < pC.size(); i++)
    {
        qxC += TMath::Cos(2 * pC[i].phi);
        qyC += TMath::Sin(2 * pC[i].phi);
    }

    qxC = qxC / pC.size();
    qyC = qyC / pC.size();
    psiC = TMath::ATan2(qyC, qxC) / 2;

    //Calculate v2
    for (unsigned int i=0; i<total_particles.size(); i++)
    {
    	float v2 = TMath::Cos(2 * (total_particles[i].phi - psiA));
    	v2s->Fill(total_particles[i].pT, v2);
    }

    res_comp->Fill(1.0, TMath::Cos(2 * (psiA - psiB)));
    res_comp->Fill(2.0, TMath::Cos(2 * (psiA - psiC)));
    res_comp->Fill(3.0, TMath::Cos(2 * (psiB - psiC)));
}

void parseampt()
{
    //Read in data file
    ifstream dataFile;
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

            //cout << partid << "    " << pv[0] << "    " << pv[1] << "    " << pv[2] << endl;

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

            //mid-rapidity
            if (p.eta > -0.5 && p.eta < 0.5)
            {
                total_particles.push_back(p);
            }

            //BBCS
            if (p.eta > -3.9 && p.eta < -3.1) pBBCS.push_back(p);

            //FVTXS
            if (p.eta > -3.1 && p.eta < -1.0) pFVTXS.push_back(p);

            //FVTXN
            if (p.eta > 1.0 && p.eta <  3.1) pFVTXN.push_back(p);

        }

        processEvent(pBBCS, pFVTXS, pFVTXN);

        pBBCS.clear();
        pFVTXS.clear();
        pFVTXN.clear();

        total_particles.clear();

        if (!dataFile) break;
    }
}

void event_plane()
{
	v2s = new TProfile("v2s", "v2s", 9, 0.2, 2, -1.0, 1.0);
    res_comp = new TProfile("res_comp", "res_comp", 3, 0.5, 3.5, -1.0, 1.0);

    //Make a file to store outputs
    TFile *fout = new TFile("eplane.root", "RECREATE");

    parseampt();

    //Save graphs
    v2s->Write();
    res_comp->Write();

    fout->Close();
}
















































