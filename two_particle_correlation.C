//-----------------------------------------------
//Code to run AMPT
//with 2 particle correlation.
//
//First:  All pT,  eta:[-3.9.-3.1] (BBCS)
//Second: pT bins, eta:[-0.5, 0.5] (Central)
//
//
//Author: P. Yin
//Date  : 08/11/16
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
vector<particle> BBCS_particle;
vector<particle> mid_rapidity_particle[9];

//-----------------------------------------------------------------------------------
//Vector declarations
int total_event_counter_amptdat = 0;

//-----------------------------------------------------------------------------------
//Graphs declarations
TH1F* event_count;

vector<TH1F*> dhis;

//-----------------------------------------------------------------------------------
//Output file


//-----------------------------------------------------------------------------------
//Test cout turn on/off set up
bool test = false;

//-----------------------------------------------------------------------------------
//Functions declarations
void processEvent()
{
	if (BBCS_particle.size() == 0) return;

	event_count->Fill(0);

	for (unsigned int first = 0; first < BBCS_particle.size(); first++)
	{
		for (int pt_bin = 0; pt_bin < 9; pt_bin++)
		{
			for (unsigned second = 0; second < mid_rapidity_particle[pt_bin].size(); second++)
			{
				float dphi = BBCS_particle[first].phi - mid_rapidity_particle[pt_bin][second].phi;

				if (dphi > 1.5 * TMath::Pi())
				{
					dphi = dphi - 2 * TMath::Pi();
				}

				if (dphi < -0.5 * TMath::Pi())
				{
					dphi = dphi + 2 * TMath::Pi();
				}

				dhis[pt_bin]->Fill(dphi);
			}
		}

		for (unsigned int second = 0; second < BBCS_particle.size(); second++)
		{
			//Skip repeated pairs
			if (first >= second) continue;

			float dphi = BBCS_particle[first].phi - BBCS_particle[second].phi;

			if (dphi > 1.5 * TMath::Pi())
			{
				dphi = dphi - 2 * TMath::Pi();
			}

			if (dphi < -0.5 * TMath::Pi())
			{
				dphi = dphi + 2 * TMath::Pi();
			}

			dhis[9]->Fill(dphi);
		}
	}
}

void parse_ampt_dat(int file_n)
{
	//Read in data file
	ifstream ampt_datFile;
	ampt_datFile.open("/Users/Bill/Desktop/Lab/TEST_Data/ampt_dau_200GeV_2fm_on_10K.dat");
	// ampt_datFile.open("/Users/Bill/Desktop/Lab/TEST_DATA/ampt_dau_200GeV_15fm_on_10K.dat");
	// ampt_datFile.open("ana/ampt.dat");

	//Skip the job if not ampt_datFile
	if (!ampt_datFile)
	{
		cout << Form("--> File ampt_%i.dat does not exist\n", file_n) << endl << endl;
		return;
	}
	else
	{
		cout << Form("--> Successfully opened file ampt_%i.dat", file_n) << endl << endl;
	}

	//In this while loop, program will read the data file line by line
	while (ampt_datFile)
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

		int Nch_counter_no_selection = 0;

		//Get the header of each event
		ampt_datFile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk;

		if (!ampt_datFile) break;

		//Analysis each particle in the event
		for (int i = 0; i < nlist; i++)
		{
			int    partid;
			float  pv[3];
			float  mass;
			double space[4];

			ampt_datFile >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

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
			p.m   = mass;
			p.x   = space[0];
			p.y   = space[1];
			p.z   = space[2];
			p.t   = space[3];

			//Store particles into vectors
			if (p.eta >= -3.9 && p.eta <= -3.1) BBCS_particle.push_back(p);
			if (p.eta >= -0.5 && p.eta <=  0.5)
			{
				for (int pt_bin = 0; pt_bin < 9; pt_bin++)
				{
					if (p.pT > 0.2 + 0.2 * pt_bin && p.pT < 0.4 + 0.2 * pt_bin)  mid_rapidity_particle[pt_bin].push_back(p);
				}
			}
		}
		processEvent();

		BBCS_particle.clear();
		for (int pt_bin = 0; pt_bin < 9; pt_bin++)
		{
			mid_rapidity_particle[pt_bin].clear();
		}
	}
}

void two_particle_correlation()
{
	event_count = new TH1F("event_count", "event_count", 1, 0, 1);

	for (int i = 0; i < 10; i++)
	{
		dhis.push_back(new TH1F(Form("d_%i", i),  "dhis", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi()));
	}

	parse_ampt_dat(0);

	TFile *fout = new TFile("two_part.root", "RECREATE");

	//Save graphs
    for(int i=0; i<10; i++)
    {   
        dhis[i]->Write();
    }

    event_count->Write();
    fout->Close();
}



















































