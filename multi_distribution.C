//-----------------------------------------------
//This code will parse AMPT final state particles
//and give distribution of:
//1. eta            distribution
//2. N_charge       distribution
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
vector<particle> projectile_nucleon_initial;
vector<particle> target_nucleon_initial;
vector<particle> parton_before_scattering;
vector<particle> nucleon_final_state;

vector<float> epsilon2;
vector<float> epsilon3;
vector<float> epsilon2_sqr;
vector<float> epsilon3_sqr;
vector<float> psi2_parton;
vector<float> v2_temp;

vector<particle> empty_vector;

//-----------------------------------------------------------------------------------
//Vector declarations
int total_event_counter_amptdat = 0;
int total_event_counter_parton_afterPropagation = 0;
int npart_sum_proj = 0;
int npart_sum_targ = 0;

//-----------------------------------------------------------------------------------
//Graphs declarations
TH1F* Nch_distribution_no_selection;
TH1F* eta_distribution_no_selection;

TH1F* event_count;

TNtuple* quark_info;
TNtuple* final_state_info;

TH1F *epsilon2_dis;
TH1F *epsilon2_sqr_dis;
TH1F *epsilon3_dis;
TH1F *epsilon3_sqr_dis;
TH1F *psi2_dis;

TH1F *v2_ebe;
TH2F *v2_ebe_epsilon2;

//-----------------------------------------------------------------------------------
//Output file
bool epsilon2_file_flag = true;
ofstream epsilon2_file;

//-----------------------------------------------------------------------------------
//Test cout turn on/off set up
bool test = false;

//-----------------------------------------------------------------------------------
//Functions declarations
// void parse_npartxy_dat(int file_n)
// {
// 	ifstream npartxy_datFile;
// 	npartxy_datFile.open("ana/npart-xy.dat");

// 	if (!npartxy_datFile)
// 	{
// 		cout << Form("--> File npart-xy_%i.dat does not exist\n", file_n) << endl << endl;
//         return;
// 	}
// 	else cout << Form("--> Successfully opened file npart-xy_%i.dat", file_n) << endl << endl;

// 	while (npartxy_datFile)
// 	{
// 		int    evtnumber;
//         int    flag;
//         int    Z_proj;
//         int    Z_targ;
//         double impactpar;

//         npartxy_datFile >> evtnumber >> flag >> Z_proj >> Z_targ >> impactpar;

//         if (!npartxy_datFile) break;

//         int total_part = Z_targ + Z_proj;

//         //Analyze each event
//         for (int i=0; i<total_part; i++)
//         {
//         	double space[3];
//         	int    sequence_number;
//         	int    status;
//         	int    present_flavor;
//         	int    original_flavor;

//         	npartxy_datFile >> space[0] >> space[1] >> sequence_number >> status >> space[2] >> present_flavor >> original_flavor;

//         	particle p;
//         	p.x = space[0];
//         	p.y = space[1];

//         	//Projectile
//         	if (sequence_number > 0)
//         	{
//         		//0: spectator;
//         		if (status == 0) continue;

//         		//1 or 2 Due to elastic collisions. 3: Due to inelastic collisions
//         		projectile_nucleon_initial.push_back(p);
//         	}

//         	//Target
//         	if (sequence_number < 0)
//         	{
//         		//0: spectator;
//         		if (status == 0) continue;

//         		//1 or 2 Due to elastic collisions. 3: Due to inelastic collisions
//         		target_nucleon_initial.push_back(p);
//         	}
//         }
// 		projectile_nucleon_initial.clear();
// 		target_nucleon_initial.clear();
// 	}
// }


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

	//Calculate epsilon_n and psi_n and put them into vectors
	//Using formulas from arXiv:1501.06880
	float e2 = TMath::Sqrt(qx2 * qx2 + qy2 * qy2) / aver2;
	float e3 = TMath::Sqrt(qx3 * qx3 + qy3 * qy3) / aver2;
	float psi2_temp = (TMath::ATan2(qy2, qx2) + TMath::Pi()) / 2;

	epsilon2.push_back(e2);
	epsilon3.push_back(e3);

	epsilon2_sqr.push_back(e2 * e2);
	epsilon3_sqr.push_back(e3 * e3);

	psi2_parton.push_back(psi2_temp);

	epsilon2_dis->Fill(e2);
	epsilon2_sqr_dis->Fill(e2 * e2);
	epsilon3_dis->Fill(e3);
	epsilon3_sqr_dis->Fill(e3 * e3);

	if (psi2_temp > 1.5 * TMath::Pi())
	{
		psi2_temp = psi2_temp - 2 * TMath::Pi();
	}

	if (psi2_temp < -0.5 * TMath::Pi())
	{
		psi2_temp = psi2_temp + 2 * TMath::Pi();
	}

	psi2_dis->Fill(psi2_temp);

	// cout << setw(10) << qx2 << setw(10) << qy2 << setw(10) << aver2 << setw(10) << e2 << setw(10) << psi2_parton << endl;
}


void parse_parton_initial_afterPropagation_dat(int file_n)
{
	//Read in data file
	ifstream partonFile;
	partonFile.open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf3/pengqi/TEST_DATA/AMPT_new_version_data/pPb/parton_afterPropagation/parton-initial-afterPropagation_%i.dat", file_n));
	// partonFile.open("ana/parton-initial-afterPropagation.dat");

	//Skip the job if not partonFile
	if (!partonFile)
	{
		cout << Form("--> File parton-initial-afterPropagation_%i.dat does not exist\n", file_n) << endl << endl;
		return;
	}
	else cout << Form("--> Successfully opened file parton-initial-afterPropagation_%i.dat", file_n) << endl << endl;

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

		if (test) cout << evtnumber << "   " << iteration << "   " << nlist << "   " << n_baryon_formed << "   " << n_meson_formed << "   " << n_inipart << "   " << n_inipart_notinzpc << endl;

		//Analysis each particle in the event
		for (int i = 0; i < nlist; i++)
		{
			int    partid;
			float  pv[3];
			float  mass;
			double space[4];

			partonFile >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

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

			if (abs(p.eta) > 3) continue;
			if (p.t > 1) continue;

			if (total_event_counter_parton_afterPropagation < 100) quark_info->Fill(p.px, p.py, p.pz, p.m, p.x, p.y, p.z, p.t, p.eta, p.pT, p.phi);

			parton_before_scattering.push_back(p);
		}
		total_event_counter_parton_afterPropagation++;

		processEvent_epsilon_psi(false, false, parton_before_scattering, empty_vector);

		parton_before_scattering.clear();
	}
}

void processEvent_v2_PP(int index)
{
	for (unsigned int i = 0; i < nucleon_final_state.size(); i++)
	{
		float v2_particle = TMath::Cos(2 * (nucleon_final_state[i].phi - psi2_parton[index]));
		v2_temp.push_back(v2_particle);
	}
}


void parse_ampt_dat(int file_n)
{
	//Read in data file
	ifstream ampt_datFile;
	ampt_datFile.open(Form("/gpfs/mnt/gpfs02/phenix/plhf/plhf3/pengqi/TEST_DATA/AMPT_new_version_data/pPb/ampt_dat/ampt_%i.dat", file_n));
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

			Nch_counter_no_selection++;

			if (total_event_counter_amptdat < 100) final_state_info->Fill(p.px, p.py, p.pz, p.m, p.x, p.y, p.z, p.t, p.eta, p.pT, p.phi);

			eta_distribution_no_selection->Fill(p.eta);

			if (abs(p.eta) < 3) nucleon_final_state.push_back(p);
		}
		processEvent_v2_PP(evtnumber - 1);

		//-----------------------------------------------
		//Calculate average v2 event by event
		float ave_v2_ebe = 0;
		for (unsigned int i=0; i<v2_temp.size(); i++)
		{
			ave_v2_ebe += v2_temp[i];
		}
		ave_v2_ebe = ave_v2_ebe / (float)v2_temp.size();
		v2_ebe->Fill(ave_v2_ebe);
		v2_ebe_epsilon2->Fill(epsilon2[total_event_counter_amptdat], ave_v2_ebe);
		//-----------------------------------------------

		nucleon_final_state.clear();
		v2_temp.clear();

		total_event_counter_amptdat++;
		event_count->Fill(0);

		npart_sum_proj += npartproj;
		npart_sum_targ += nparttarg;

		Nch_distribution_no_selection->Fill(Nch_counter_no_selection);
	}
}

void multi_distribution()
{
	event_count = new TH1F("event_count", "event_count", 1, 0, 1);

	Nch_distribution_no_selection = new TH1F("Nch_distribution_no_selection", "Nch_distribution_no_selection", 500, -0.5, 999.5);
	eta_distribution_no_selection = new TH1F("eta_distribution_no_selection", "eta_distribution_no_selection", 200, -10, 10);

	epsilon2_dis = new TH1F("epsilon2_dis", "epsilon2_dis", 100, -0.01, 0.99);
	epsilon2_sqr_dis = new TH1F("epsilon2_sqr_dis", "epsilon2_sqr_dis", 100, -0.01, 0.99);
	epsilon3_dis = new TH1F("epsilon3_dis", "epsilon3_dis", 100, -0.01, 0.99);
	epsilon3_sqr_dis = new TH1F("epsilon3_sqr_dis", "epsilon3_sqr_dis", 100, -0.01, 0.99);
	psi2_dis = new TH1F("psi2_dis", "psi2_dis", 50, -0.5 * TMath::Pi(), 1.5 * TMath::Pi());
	v2_ebe = new TH1F("v2_ebe", "v2_ebe", 100, 0, 1);
	v2_ebe_epsilon2 = new TH2F("v2_ebe_epsilon2", "v2_ebe_epsilon2", 200, -0.01, 0.99, 100, -0.01, 0.49);

	quark_info = new TNtuple("quark_info", "quark_info", "px:py:pz:mass:x:y:z:t:eta:pT:phi", 32000);
	final_state_info = new TNtuple("final_state_info", "final_state_info", "px:py:pz:mass:x:y:z:t:eta:pT:phi", 32000);

	for (int i = 0; i < 100; i++)
	{
		parse_parton_initial_afterPropagation_dat(i);
		parse_ampt_dat(i);
		psi2_parton.clear();
		// parse_npartxy_dat(i);
	}

	float ep2avg = 0;
	float ep3avg = 0;
	float ep2_sqr_avg = 0;
	float ep3_sqr_avg = 0;

	for (unsigned int i = 0; i < epsilon2.size(); i++)
	{
		ep2avg += epsilon2[i];
	}
	ep2avg = ep2avg / epsilon2.size();

	for (unsigned int i = 0; i < epsilon3.size(); i++)
	{
		ep3avg += epsilon3[i];
	}
	ep3avg = ep3avg / epsilon3.size();

	for (unsigned int i = 0; i < epsilon2_sqr.size(); i++)
	{
		ep2_sqr_avg += epsilon2_sqr[i];
	}
	ep2_sqr_avg = ep2_sqr_avg / epsilon2_sqr.size();

	for (unsigned int i = 0; i < epsilon3_sqr.size(); i++)
	{
		ep3_sqr_avg += epsilon3_sqr[i];
	}
	ep3_sqr_avg = ep3_sqr_avg / epsilon3_sqr.size();

	Nch_distribution_no_selection->Scale(1. / ((float)total_event_counter_amptdat * Nch_distribution_no_selection->GetBinWidth(1)));
	eta_distribution_no_selection->Scale(1. / ((float)total_event_counter_amptdat * eta_distribution_no_selection->GetBinWidth(1)));

	TFile *fout = new TFile("results.root", "RECREATE");

	event_count->Write();

	Nch_distribution_no_selection->Write();
	eta_distribution_no_selection->Write();

	epsilon2_dis->Write();
	epsilon2_sqr_dis->Write();
	epsilon3_dis->Write();
	epsilon3_sqr_dis->Write();
	psi2_dis->Write();
	v2_ebe->Write();
	v2_ebe_epsilon2->Write();

	quark_info->Write();
	final_state_info->Write();

	fout->Close();

	cout << "------------- Output -------------" << endl;
	cout << "<epsilon2>         | " << ep2avg << endl;
	cout << "----------------------------------" << endl;
	cout << "<epsilon2 square>  | " << ep2_sqr_avg << endl;
	cout << "----------------------------------" << endl;
	cout << "<epsilon3>         | " << ep3avg << endl;
	cout << "----------------------------------" << endl;
	cout << "<epsilon3 square>  | " << ep3_sqr_avg << endl;
	cout << "----------------------------------" << endl;
	cout << "<Npart_targ>       | " << (float)npart_sum_targ / (float)total_event_counter_amptdat << endl;
	cout << "----------------------------------" << endl;
	cout << "<Npart_proj>       | " << (float)npart_sum_proj / (float)total_event_counter_amptdat << endl;
	cout << "----------------------------------" << endl;
	cout << "<Npart>            | " << ((float)npart_sum_proj + (float)npart_sum_targ) / (float)total_event_counter_amptdat << endl;
	cout << "----------------------------------" << endl << endl;

	if (epsilon2_file_flag)
	{
		epsilon2_file.open("epsilon2.txt");
		for (unsigned int i=0; i<epsilon2.size(); i++)
		{
			epsilon2_file << epsilon2[i] << endl;
		}
		epsilon2_file.close();
	}
}







































