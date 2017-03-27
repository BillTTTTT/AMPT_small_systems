//-----------------------------------------------
// Event Plane Method
//
// 3-sub event method resolution from:
// BBC, FVTXS, CNT particles.
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

vector<particle> pA;
vector<particle> pB;
vector<particle> pC;

TProfile* v2s_FVTXS[6];
TProfile* res_comp[6];

TH1F* dhis_bbcs;
TH1F* hcount;

TH2D *eff_fvtx_s;
TH2D *eff_fvtx_n;

TF1 *frandom;

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
    for (unsigned int i=0; i<pC.size(); i++)
    {
        float v2 = TMath::Cos(2 * (pC[i].phi - psiA));
        v2s_FVTXS[index]->Fill(pC[i].pT, v2);
    }

    res_comp[index]->Fill(1.0, TMath::Cos(2 * (psiA - psiB)));
    res_comp[index]->Fill(2.0, TMath::Cos(2 * (psiA - psiC)));
    res_comp[index]->Fill(3.0, TMath::Cos(2 * (psiB - psiC)));
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

            //mid-rapidity
            if (ifCNT(p.eta))
            {
                pC.push_back(p);
            }

            //BBCS
            if (ifBBCS(p.eta)) 
            {
                pB.push_back(p);
                ct_bbcs++;
                hcount->Fill(0);
            }

            //FVTXS
            if (ifFVTXS(p.eta) && test_eff_s(p.pT, p.eta)) 
            {
                pA.push_back(p);
                ct_fvtxs++;
                hcount->Fill(1);
            }
        }

        dhis_bbcs->Fill(ct_bbcs);

        if (ct_bbcs >= 28) processEvent(0);
        if (ct_bbcs >= 24 && ct_bbcs < 28) processEvent(1);
        if (ct_bbcs >= 19 && ct_bbcs < 24) processEvent(2);
        if (ct_bbcs >= 12 && ct_bbcs < 19) processEvent(3);
        if (ct_bbcs >=  7 && ct_bbcs < 12) processEvent(4);
        if (ct_bbcs < 7) processEvent(5);

        pA.clear();
        pB.clear();
        pC.clear();

        hcount->Fill(2);

        if (!dataFile) break;
    }
}

void event_plane()
{
    const int NPTBINS = 16;
    double ptlim[NPTBINS + 1] = {
        0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
        2.5, 3.0, 3.5, 4.0, 4.5, 5.0
    };
    
    v2s_FVTXS[0] = new TProfile("v2s_FVTXS_0", "v2s_FVTXS_0", NPTBINS, ptlim, -1.1, 1.1);
    res_comp[0] = new TProfile("res_comp_0", "res_comp_0", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_FVTXS[1] = new TProfile("v2s_FVTXS_1", "v2s_FVTXS_1", NPTBINS, ptlim, -1.1, 1.1);
    res_comp[1] = new TProfile("res_comp_1", "res_comp_1", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_FVTXS[2] = new TProfile("v2s_FVTXS_2", "v2s_FVTXS_2", NPTBINS, ptlim, -1.1, 1.1);
    res_comp[2] = new TProfile("res_comp_2", "res_comp_2", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_FVTXS[3] = new TProfile("v2s_FVTXS_3", "v2s_FVTXS_3", NPTBINS, ptlim, -1.1, 1.1);
    res_comp[3] = new TProfile("res_comp_3", "res_comp_3", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_FVTXS[4] = new TProfile("v2s_FVTXS_4", "v2s_FVTXS_4", NPTBINS, ptlim, -1.1, 1.1);
    res_comp[4] = new TProfile("res_comp_4", "res_comp_4", 6, 0.5, 6.5, -1.0, 1.0);

    v2s_FVTXS[5] = new TProfile("v2s_FVTXS_5", "v2s_FVTXS_5", NPTBINS, ptlim, -1.1, 1.1);
    res_comp[5] = new TProfile("res_comp_5", "res_comp_5", 6, 0.5, 6.5, -1.0, 1.0);

    dhis_bbcs = new TH1F("dhis_bbcs", "dhis_bbcs", 200, -0.5, 199.5);

    hcount = new TH1F("hcount", "count", 6, 0, 6);

    TFile *f_fvtxs = new TFile("fvtx_acc.root");
    TFile *f_fvtxn = new TFile("fvtx_acc_n.root");

    eff_fvtx_s = (TH2D*)f_fvtxs->Get("rh");
    eff_fvtx_n = (TH2D*)f_fvtxn->Get("rh");

    frandom = new TF1("frandom", "1", 0.0, 1.0);

    //Make a file to store outputs
    TFile *fout = new TFile("out_EP.root", "RECREATE");

    parseampt();

    //Save graphs
    for (int i=0; i<6; i++)
    {
        v2s_FVTXS[i]->Write();
        res_comp[i]->Write();
    }
    hcount->Write();
    dhis_bbcs->Write();
    
    fout->Close();
}


























