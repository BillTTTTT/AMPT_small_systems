//-----------------------------------------------
//Code to run AMPT d+Au @ 200GeV
//with event plane method
//for particles in BBC, FVTXS, FVTXN, particles.
//using 3-sub event method for resolution calc.
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

//-----------------------------------------------
//Variables
//-----------------------------------------------

struct particle
{
    int   id;
    float px;
    float py;
    float pz;
    float x;
    float y;
    float z;
    float eta;
    float phi;
    float pT;
};

//Number of nucleons in system
// --> p+Au = 198
// --> d+Au = 199
// --> d+Pb = 210
// --> p+Pb = 209

//BBC [-3.9,-3)
vector<particle> pA;
vector<particle> total_particles;

//FVTXS [-3,-1)
vector<particle> pB;

//FVTXN [1,3)
vector<particle> pC;

TProfile* v2s_BBCS[6]; //v2s_BBCS
// TProfile* v2s_FVTXS; //v2s_BBCS
TProfile* res_comp[6];
// TH1F* dhis_1;
// TH1F* dhis_2;
// TH1F* dhis_3;

TH1F* dhis_bbcs;
TH1F* dhis_fvtxs;

TH1F* eta_distribution_no_selection;

TH1F *hcount;

void  processEvent(int index)
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
    // cout << psiA << endl;

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
        // float v2_fvtxs = TMath::Cos(2 * (total_particles[i].phi - psiB));
        // if (total_particles[i].pT>0.2 && total_particles[i].pT<0.4) cout << v2 << endl;
        v2s_BBCS[index]->Fill(total_particles[i].pT, v2);
        // v2s_FVTXS->Fill(total_particles[i].pT, v2_fvtxs);
    }

    res_comp[index]->Fill(1.0, TMath::Cos(2 * (psiA - psiB)));
    res_comp[index]->Fill(2.0, TMath::Cos(2 * (psiA - psiC)));
    res_comp[index]->Fill(3.0, TMath::Cos(2 * (psiB - psiC)));

    // res_comp->Fill(4.0, TMath::Cos(2 * (psiB - psiA)));
    // res_comp->Fill(5.0, TMath::Cos(2 * (psiB - psiC)));
    // res_comp->Fill(6.0, TMath::Cos(2 * (psiA - psiC)));

    // dhis_1->Fill(psiA - psiB);
    // dhis_2->Fill(psiA - psiC);
    // dhis_3->Fill(psiB - psiC);
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

        int ct_bbcs = 0;
        int ct_fvtxs = 0;

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

            eta_distribution_no_selection->Fill(p.eta);

            //mid-rapidity
            if (p.eta > -0.35 && p.eta < 0.35)
            {
                total_particles.push_back(p);
            }

            //BBCS
            if (p.eta > -3.9 && p.eta < -3.1) 
            {
                pA.push_back(p);
                ct_bbcs++;
                hcount->Fill(0);
            }

            //FVTXS
            if (p.eta > -3.1 && p.eta < -1.0) 
            {
                pB.push_back(p);
                ct_fvtxs++;
                hcount->Fill(1);
            }

            //FVTXN
            if (p.eta > 1.0 && p.eta <  3.1) pC.push_back(p);

        }

        dhis_bbcs->Fill(ct_bbcs);
        // dhis_fvtxs->Fill(ct_fvtxs);

        // cout << total_particles.size() << endl;
        if (ct_bbcs >= 28) processEvent(0);
        if (ct_bbcs >= 24 && ct_bbcs < 28) processEvent(1);
        if (ct_bbcs >= 19 && ct_bbcs < 24) processEvent(2);
        if (ct_bbcs >= 12 && ct_bbcs < 19) processEvent(3);
        if (ct_bbcs >=  7 && ct_bbcs < 12) processEvent(4);
        if (ct_bbcs < 7) processEvent(5);

        pA.clear();
        pB.clear();
        pC.clear();

        total_particles.clear();

        hcount->Fill(2);

        if (!dataFile) break;
    }
}

void f1_ep_BBCS_0()
{
    v2s_BBCS[0] = new TProfile("v2s_BBCS_0", "v2s_BBCS_0", 14, 0.2, 3, -1.0, 1.0);
    res_comp[0] = new TProfile("res_comp_0", "res_comp_0", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_BBCS[1] = new TProfile("v2s_BBCS_1", "v2s_BBCS_1", 14, 0.2, 3, -1.0, 1.0);
    res_comp[1] = new TProfile("res_comp_1", "res_comp_1", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_BBCS[2] = new TProfile("v2s_BBCS_2", "v2s_BBCS_2", 14, 0.2, 3, -1.0, 1.0);
    res_comp[2] = new TProfile("res_comp_2", "res_comp_2", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_BBCS[3] = new TProfile("v2s_BBCS_3", "v2s_BBCS_3", 14, 0.2, 3, -1.0, 1.0);
    res_comp[3] = new TProfile("res_comp_3", "res_comp_3", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_BBCS[4] = new TProfile("v2s_BBCS_4", "v2s_BBCS_4", 14, 0.2, 3, -1.0, 1.0);
    res_comp[4] = new TProfile("res_comp_4", "res_comp_4", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_BBCS[5] = new TProfile("v2s_BBCS_5", "v2s_BBCS_5", 14, 0.2, 3, -1.0, 1.0);
    res_comp[5] = new TProfile("res_comp_5", "res_comp_5", 6, 0.5, 6.5, -1.0, 1.0);
    // v2s_FVTXS = new TProfile("v2s_FVTXS", "v2s_FVTXS", 14, 0.2, 3, -1.0, 1.0);
    // dhis_1 = new TH1F("dhis_1", "dhis_1", 50, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    // dhis_2 = new TH1F("dhis_2", "dhis_2", 50, -0.5*TMath::Pi(), 1.5*TMath::Pi());
    // dhis_3 = new TH1F("dhis_3", "dhis_3", 50, -0.5*TMath::Pi(), 1.5*TMath::Pi());

    dhis_bbcs = new TH1F("dhis_bbcs", "dhis_bbcs", 200, -0.5, 199.5);
    // dhis_fvtxs = new TH1F("dhis_fvtxs", "dhis_fvtxs", 200, -0.5, 199.5);

    eta_distribution_no_selection = new TH1F("eta_distribution_no_selection", "eta_distribution_no_selection", 200, -10, 10);

    hcount = new TH1F("hcount", "count", 6, 0, 6);

    //Make a file to store outputs
    TFile *fout = new TFile("ep_bbcs_0.root", "RECREATE");

    parseampt();

    //Save graphs

    for (int i=0; i<6; i++)
    {
        v2s_BBCS[i]->Write();
        res_comp[i]->Write();
    }
    
    // v2s_FVTXS->Write();
    // dhis_1->Write();
    // dhis_2->Write();
    // dhis_3->Write();
    hcount->Write();

    eta_distribution_no_selection->Write();

    dhis_bbcs->Write();
    // dhis_fvtxs->Write();
    
    fout->Close();
}


























