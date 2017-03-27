//-----------------------------------------------
// Parton Participant Plane Method
//
// Author: P. Yin
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

#include "range.h"

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
vector<float> ep2;

vector<particle> final_p;

//-----------------------------------------------------------------------------------
//Graphs declarations
TProfile* v2s;
TProfile* v2s_pt[6];

TProfile* epsilon2_nch;

TH1F *epsilon2_dis;
TH1F *dhis_v2;
TH1F *dhis_qn;
TH1F *dhis_b;
TH1F *dhis_eta;

TH2D *eff_fvtx_s;
TH2D *eff_fvtx_n;

TF1 *frandom;

//-----------------------------------------------------------------------------------
//Test cout turn on/off set up
bool test = false;

//-----------------------------------------------------------------------------------
//Functions declarations
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

	// qx3 = qx3 / (float)count;
	// qy3 = qy3 / (float)count;

	aver2    = aver2 / (float)count;

	float temp_psi2 = (TMath::ATan2(qy2, qx2) + TMath::Pi()) / 2;
	float e2 = TMath::Sqrt(qx2 * qx2 + qy2 * qy2) / aver2;

	epsilon2_dis->Fill(e2);

	ep2.push_back(e2);
	psi2.push_back(temp_psi2);
}


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


void  processEvent_ampt(int evtnumber, int ncharge, int index)
{
	//Calculate v2
	for (unsigned int i=0; i<final_p.size(); i++)
	{
		float v2 = TMath::Cos(2 * (final_p[i].phi - psi2[evtnumber]));
		v2s->Fill(ncharge, v2);
		v2s_pt[index]->Fill(final_p[i].pT, v2);
		
		dhis_v2->Fill(v2);
	}

	epsilon2_nch->Fill(ncharge, ep2[evtnumber]);

	float cmx = 0;
	float cmy = 0;
	float qx2 = 0;
	float qy2 = 0;
	float aver2 = 0;

	for (unsigned int i = 0; i < final_p.size(); i++)
	{
		cmx = cmx + final_p[i].x;
		cmy = cmy + final_p[i].y;
	}

	int count = final_p.size();

	cmx = cmx / (float)count;
	cmy = cmy / (float)count;

	for (unsigned int i = 0; i < final_p.size(); i++)
	{
		//Shift to center of mass frame
		final_p[i].x = final_p[i].x - cmx;
		final_p[i].y = final_p[i].y - cmy;

		final_p[i].phi = TMath::ATan2(final_p[i].y, final_p[i].x);
	}

	for (unsigned int i = 0; i < final_p.size(); i++)
	{
		qx2 = qx2 + TMath::Cos(2 * final_p[i].phi);
		qy2 = qy2 + TMath::Sin(2 * final_p[i].phi);
	}

	float numerator = TMath::Sqrt(qx2 * qx2 + qy2 * qy2);
	float denominator = TMath::Sqrt(final_p.size());

	float qn = numerator / denominator;

	dhis_qn->Fill(qn);
}

bool test_eff_s(float pT, float eta)
{
    int pTbin = eff_fvtx_s->GetXaxis()->FindBin(pT);
    int etabin = eff_fvtx_s->GetYaxis()->FindBin(eta);

    float n = eff_fvtx_s->GetBinContent(pTbin, etabin);

    float test = frandom->GetRandom();

    if (test < n) return true;
    else return false;
}

bool test_eff_n(float pT, float eta)
{
    int pTbin = eff_fvtx_n->GetXaxis()->FindBin(pT);
    int etabin = eff_fvtx_n->GetYaxis()->FindBin(eta);

    float n = eff_fvtx_n->GetBinContent(pTbin, etabin);

    float test = frandom->GetRandom();

    if (test < n) return true;
    else return false;
}


void parse_ampt_file()
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

        int    ncharge = 0;
        int ct_bbcs = 0;

        //Get the header of each event
        dataFile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk;

        if (!dataFile) break;

        dhis_b->Fill(impactpar);

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

            dhis_eta->Fill(p.eta);

            if (ifCNT(p.eta)) 
            {
            	final_p.push_back(p);
            	ncharge++;
            }

            if (ifBBCS(p.eta)) ct_bbcs++;

        }

        if (ct_bbcs >= 28) processEvent_ampt(event_counter, ncharge, 0);
        if (ct_bbcs >= 24 && ct_bbcs < 28) processEvent_ampt(event_counter, ncharge, 1);
        if (ct_bbcs >= 19 && ct_bbcs < 24) processEvent_ampt(event_counter, ncharge, 2);
        if (ct_bbcs >= 12 && ct_bbcs < 19) processEvent_ampt(event_counter, ncharge, 3);
        if (ct_bbcs >=  7 && ct_bbcs < 12) processEvent_ampt(event_counter, ncharge, 4);
        if (ct_bbcs < 7) processEvent_ampt(event_counter, ncharge, 5);

        final_p.clear();

        event_counter += 1;

        if (!dataFile) break;
    }
}


void parton_pplane()
{
    TFile *fout = new TFile("ppplane.root", "RECREATE");
  
    const int NPTBINS = 16;
	double ptlim[NPTBINS + 1] = {
		0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
		2.5, 3.0, 3.5, 4.0, 4.5, 5.0
	};

    // 100, -0.5, 5999.5, -10, 10 Pb+Pb
    // 50, -0.5, 199.5, -10, 10    d+Au
    v2s = new TProfile("v2s", "v2s", 50, -0.5, 199.5, -10, 10); //
    v2s_pt[0] = new TProfile("v2s_pt_0", "v2s_pt_0", NPTBINS, ptlim, -1.1, 1.1);
    v2s_pt[1] = new TProfile("v2s_pt_1", "v2s_pt_1", NPTBINS, ptlim, -1.1, 1.1);
    v2s_pt[2] = new TProfile("v2s_pt_2", "v2s_pt_2", NPTBINS, ptlim, -1.1, 1.1);
    v2s_pt[3] = new TProfile("v2s_pt_3", "v2s_pt_3", NPTBINS, ptlim, -1.1, 1.1);
    v2s_pt[4] = new TProfile("v2s_pt_4", "v2s_pt_4", NPTBINS, ptlim, -1.1, 1.1);
    v2s_pt[5] = new TProfile("v2s_pt_5", "v2s_pt_5", NPTBINS, ptlim, -1.1, 1.1);

    epsilon2_nch = new TProfile("epsilon2_nch", "epsilon2_nch", 100, -0.5, 1499.5, -10, 10); //

    dhis_v2 = new TH1F("dhis_v2", "dhis_v2", 200, -1, 1);
    epsilon2_dis = new TH1F("epsilon2_dis", "epsilon2_dis", 100, -0.01, 0.99);
    dhis_qn = new TH1F("dhis_qn", "dhis_qn", 200, 0, 20);
    dhis_b = new TH1F("dhis_b", "dhis_b", 50, 0, 20);
    dhis_eta = new TH1F("dhis_eta", "dhis_eta", 100, -10, 10);

    TFile *f_fvtxs = new TFile("fvtx_acc.root");
    TFile *f_fvtxn = new TFile("fvtx_acc_n.root");

    eff_fvtx_s = (TH2D*)f_fvtxs->Get("rh");
    eff_fvtx_n = (TH2D*)f_fvtxn->Get("rh");

    frandom = new TF1("frandom", "1", 0.0, 1.0);

    parse_afterPropagation_file();
    parse_ampt_file();
    
    fout->Write();
    fout->Close();
}













































