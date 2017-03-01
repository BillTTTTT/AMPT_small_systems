//-----------------------------------------------
//Code to run AMPT
//participant plane method with nucleons
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
int total_event_counter_amptdat = 0;

//-----------------------------------------------------------------------------------
//Vector declarations
vector<particle> projectile_nucleon_initial;
vector<particle> target_nucleon_initial;
vector<particle> nucleon_final_state[9];

vector<float> psi2_nucleon;

//-----------------------------------------------------------------------------------
//Graphs declarations
TProfile* v2s;

//-----------------------------------------------------------------------------------
//Output file

//-----------------------------------------------------------------------------------
//Test cout turn on/off set up
bool test = false;

//-----------------------------------------------------------------------------------
//Functions declarations
void processEvent_epsilon_psi(bool WithProj, bool Smear, vector<particle> p_targ, vector<particle> p_proj)
{
	float cmx = 0;
	float cmy = 0;
	float qx2 = 0;
	float qx3 = 0;
	float qy2 = 0;
	float qy3 = 0;
	float aver2 = 0;

	TF1 *fgous = new TF1("fgous", "x*TMath::Exp(-(x*x)/(2*[0]*[0]))", 0.0, 2.0);
	TF1 *fphi  = new TF1("fphi", "1.0", 0.0, 2.0 * TMath::Pi());
	fgous->FixParameter(0, 0.4);

	//Calculate centroid
	for (unsigned int i = 0; i < p_targ.size(); i++)
	{
		cmx = cmx + p_targ[i].x;
		cmy = cmy + p_targ[i].y;
	}

	int count = p_targ.size();

	if (WithProj)
	{
		for (unsigned int i = 0; i < p_proj.size(); i++)
		{
			cmx = cmx + p_proj[i].x;
			cmy = cmy + p_proj[i].y;
		}
		count = count + p_proj.size();
	}

	cmx = cmx / (float)count;
	cmy = cmy / (float)count;

	// cout << cmx << "   " << cmy << "   " << count << endl;

	//Calculate each useful value of collision particles
	//Store them in vectors
	for (unsigned int i = 0; i < p_targ.size(); i++)
	{
		//Shift to center of mass frame
		p_targ[i].x = p_targ[i].x - cmx;
		p_targ[i].y = p_targ[i].y - cmy;

		p_targ[i].phi = TMath::ATan2(p_targ[i].y, p_targ[i].x);
		p_targ[i].rsquare = p_targ[i].x * p_targ[i].x + p_targ[i].y * p_targ[i].y;
	}

	if (WithProj)
	{
		for (unsigned int i = 0; i < p_proj.size(); i++)
		{
			//Shift to center of mass frame
			p_proj[i].x = p_proj[i].x - cmx;
			p_proj[i].y = p_proj[i].y - cmy;

			p_proj[i].phi = TMath::ATan2(p_proj[i].y, p_proj[i].x);
			p_proj[i].rsquare = p_proj[i].x * p_proj[i].x + p_proj[i].y * p_proj[i].y;
		}
	}

	if (Smear) count = 0;

	//Calculate the average values for computing epsilon_2 & 3
	for (unsigned int i = 0; i < p_targ.size(); i++)
	{
		if (!Smear)
		{
			qx2 = qx2 + p_targ[i].rsquare * TMath::Cos(2 * p_targ[i].phi);
			qy2 = qy2 + p_targ[i].rsquare * TMath::Sin(2 * p_targ[i].phi);

			qx3 = qx3 + p_targ[i].rsquare * TMath::Cos(3 * p_targ[i].phi);
			qy3 = qy3 + p_targ[i].rsquare * TMath::Sin(3 * p_targ[i].phi);

			aver2    = aver2 + p_targ[i].rsquare;
		}

		if (Smear)
		{
			for (int j = 0; j < 100; j++)
			{
				float random_r = fgous->GetRandom();
				float random_phi = fphi->GetRandom();

				float x_temp = p_targ[i].x + random_r * TMath::Cos(random_phi);
				float y_temp = p_targ[i].y + random_r * TMath::Sin(random_phi);

				float rsquare_temp = x_temp * x_temp + y_temp * y_temp;
				float phi_temp = TMath::ATan2(y_temp, x_temp);

				qx2 = qx2 + rsquare_temp * TMath::Cos(2 * phi_temp);
				qy2 = qy2 + rsquare_temp * TMath::Sin(2 * phi_temp);

				qx3 = qx3 + rsquare_temp * TMath::Cos(3 * phi_temp);
				qy3 = qy3 + rsquare_temp * TMath::Sin(3 * phi_temp);

				aver2    = aver2 + rsquare_temp;

				count = count + 1;
			}
		}
	}

	if (WithProj)
	{
		for (unsigned int i = 0; i < p_proj.size(); i++)
		{
			if (!Smear)
			{
				qx2 = qx2 + p_proj[i].rsquare * TMath::Cos(2 * p_proj[i].phi);
				qy2 = qy2 + p_proj[i].rsquare * TMath::Sin(2 * p_proj[i].phi);

				qx3 = qx3 + p_proj[i].rsquare * TMath::Cos(3 * p_proj[i].phi);
				qy3 = qy3 + p_proj[i].rsquare * TMath::Sin(3 * p_proj[i].phi);

				aver2    = aver2 + p_proj[i].rsquare;
			}

			if (Smear)
			{
				for (int j = 0; j < 100; j++)
				{
					float random_r = fgous->GetRandom();
					float random_phi = fphi->GetRandom();

					float x_temp = p_proj[i].x + random_r * TMath::Cos(random_phi);
					float y_temp = p_proj[i].y + random_r * TMath::Sin(random_phi);

					float rsquare_temp = x_temp * x_temp + y_temp * y_temp;
					float phi_temp = TMath::ATan2(y_temp, x_temp);

					qx2 = qx2 + rsquare_temp * TMath::Cos(2 * phi_temp);
					qy2 = qy2 + rsquare_temp * TMath::Sin(2 * phi_temp);

					qx3 = qx3 + rsquare_temp * TMath::Cos(3 * phi_temp);
					qy3 = qy3 + rsquare_temp * TMath::Sin(3 * phi_temp);

					aver2    = aver2 + rsquare_temp;

					count = count + 1;
				}
			}
		}
	}

	qx2 = qx2 / (float)count;
	qy2 = qy2 / (float)count;

	qx3 = qx3 / (float)count;
	qy3 = qy3 / (float)count;

	aver2    = aver2 / (float)count;

	float psi2_temp = (TMath::ATan2(qy2, qx2) + TMath::Pi()) / 2;

	psi2_nucleon.push_back(psi2_temp);
}

void parse_npartxy_dat(int file_n)
{
	ifstream npartxy_datFile;
	npartxy_datFile.open("npart-xy.dat");

	if (!npartxy_datFile)
	{
		cout << Form("--> File npart-xy_%i.dat does not exist\n", file_n) << endl << endl;
		return;
	}
	else cout << Form("--> Successfully opened file npart-xy_%i.dat", file_n) << endl << endl;

	while (npartxy_datFile)
	{
		int    evtnumber;
		int    flag;
		int    Z_proj;
		int    Z_targ;
		double impactpar;

		npartxy_datFile >> evtnumber >> flag >> Z_proj >> Z_targ >> impactpar;

		if (!npartxy_datFile) break;

		int total_part = Z_targ + Z_proj;

		//Analyze each event
		for (int i = 0; i < total_part; i++)
		{
			double space[3];
			int    sequence_number;
			int    status;
			int    present_flavor;
			int    original_flavor;

			npartxy_datFile >> space[0] >> space[1] >> sequence_number >> status >> space[2] >> present_flavor >> original_flavor;

			particle p;
			p.x = space[0];
			p.y = space[1];

			//Projectile
			if (sequence_number > 0)
			{
				//0: spectator;
				if (status == 0) continue;

				//1 or 2 Due to elastic collisions. 3: Due to inelastic collisions
				projectile_nucleon_initial.push_back(p);
			}

			//Target
			if (sequence_number < 0)
			{
				//0: spectator;
				if (status == 0) continue;

				//1 or 2 Due to elastic collisions. 3: Due to inelastic collisions
				target_nucleon_initial.push_back(p);
			}
		}
		processEvent_epsilon_psi(true, true, target_nucleon_initial, projectile_nucleon_initial);

		projectile_nucleon_initial.clear();
		target_nucleon_initial.clear();
	}
}

void  processEvent()
{
	//Calculate v2
	for (int i = 0; i < 9; i++)
	{
		for (unsigned int j = 0; j < nucleon_final_state[i].size(); j++)
		{
			float v2 = TMath::Cos(2 * (nucleon_final_state[i][j].phi - psi2_nucleon[total_event_counter_amptdat]));
			v2s->Fill(0.3 + i * 0.2, v2);
		}
	}
}

void parse_ampt_dat(int file_n)
{
	//Read in data file
	ifstream ampt_datFile;
	ampt_datFile.open("ampt.dat");
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
			if (p.eta >= -0.5 && p.eta <=  0.5)
			{
				for (int pt_bin = 0; pt_bin < 9; pt_bin++)
				{
					if (p.pT > 0.2 + 0.2 * pt_bin && p.pT < 0.4 + 0.2 * pt_bin)  nucleon_final_state[pt_bin].push_back(p);
				}
			}
		}
		processEvent();
		total_event_counter_amptdat++;
		for (int pt_bin = 0; pt_bin < 9; pt_bin++)
		{
			nucleon_final_state[pt_bin].clear();
		}
	}
}

void nucleon_pplane()
{
	v2s = new TProfile("v2s", "v2s", 10, 0, 2, -10, 10);

	parse_npartxy_dat(0);
	parse_ampt_dat(0);

	//Make a file to store outputs
    TFile *fout = new TFile("nucleon_participant.root", "RECREATE");

	//Save graphs
    v2s->Write();

    fout->Close();
}








































