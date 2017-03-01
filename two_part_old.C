//-----------------------------------------------
//Code to run AMPT @ 20 GeV and 200GeV
//with 2 particle correlation.
//First:  All pT,  eta:[-3.9.-3.1] (BBC)
//Second: pT bins, eta:[-0.35, 0.35] (Central)
//Also contains a 6th graph [-3.9.-3.1][-3.9.-3.1]
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

vector<particle> p1;
vector<particle> p21; //pT:[0.2,0.4)
vector<particle> p22; //pT:[0.4,0.6)
vector<particle> p23; //pT:[0.6,0.8)
vector<particle> p24; //pT:[0.8,1.0)
vector<particle> p25; //pT:[1.0,1.2)
vector<particle> p26; //pT:[1.2,1.4)
vector<particle> p27; //pT:[1.4,1.6)
vector<particle> p28; //pT:[1.6,1.8)
vector<particle> p29; //pT:[1.8,2.0)

vector<TH1F*> dhis;

TH1F *hcount;

void processEvent()
{
    if(p1.size() == 0) return;

    hcount->Fill(0);

    for(unsigned int i=0; i<p1.size(); i++)
    {
        for(unsigned int j=0; j<p21.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p21[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[0]->Fill(dphi);
        }

        for(unsigned int j=0; j<p22.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p22[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[1]->Fill(dphi);
        }

        for(unsigned int j=0; j<p23.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p23[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[2]->Fill(dphi);
        }

        for(unsigned int j=0; j<p24.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p24[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[3]->Fill(dphi);
        }

        for(unsigned int j=0; j<p25.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p25[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[4]->Fill(dphi);
        }

        for(unsigned int j=0; j<p26.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p26[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[5]->Fill(dphi);
        }

        for(unsigned int j=0; j<p27.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p27[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[6]->Fill(dphi);
        }

        for(unsigned int j=0; j<p28.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p28[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[7]->Fill(dphi);
        }

        for(unsigned int j=0; j<p29.size(); j++)
        {
            float dphi;

            dphi = p1[i].phi - p29[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[8]->Fill(dphi);
        }

        for(unsigned int j=0; j<p1.size(); j++)
        {
            //Skip repeated pairs
            if(i>=j) continue;

            float dphi;

            dphi = p1[i].phi - p1[j].phi;

            if(dphi > 1.5*TMath::Pi())
            {
                dphi = dphi - 2 * TMath::Pi();
            }

            if(dphi < -0.5*TMath::Pi())
            {
                dphi = dphi + 2 * TMath::Pi();
            }

            dhis[9]->Fill(dphi);
        }
    }
}

void parseampt()
{
    //Read in data file
    ifstream dataFile;
    dataFile.open("ana/ampt.dat");

    //Skip the job if not dataFile
    if(!dataFile)  return;

    //In this while loop, program will read the data file line by line
    while(dataFile)
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

        if(!dataFile) break;

        //Analysis each particle in the event
        for (int i=0; i<nlist; i++)
        {
            int    partid;
            float  pv[3];
            float  mass;
            double space[4];

            dataFile >> partid >> pv[0] >> pv[1] >> pv[2] >> mass >> space[0] >> space[1] >> space[2] >> space[3];

            //Skip non-charged particles that we are not interested
            //
            //+-211 are pions,  +-321 are kaons, +-2212 are protons
            if(abs(partid) != 211 && abs(partid) != 321 && abs(partid) != 2212) continue;

            if (TMath::Sqrt(pv[0] * pv[0] + pv[1] * pv[1]) < 0.0001) continue;

            //Calculate the energy
            float energy = TMath::Sqrt(pv[0]*pv[0] + pv[1]*pv[1] + pv[2]*pv[2] + mass*mass);

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
            p.x   = space[0];
            p.y   = space[1];
            p.z   = space[2];

            //Store particles into vectors
            if(p.eta >= -3.9 && p.eta <= -3.1) p1.push_back(p);
            if(p.eta >= -0.35 && p.eta <=  0.35)
            {
                if(p.pT >= 0.2 && p.pT < 0.4) p21.push_back(p);
                if(p.pT >= 0.4 && p.pT < 0.6) p22.push_back(p);
                if(p.pT >= 0.6 && p.pT < 0.8) p23.push_back(p);
                if(p.pT >= 0.8 && p.pT < 1.0) p24.push_back(p);
                if(p.pT >= 1.0 && p.pT < 1.2) p25.push_back(p);
                if(p.pT >= 1.2 && p.pT < 1.4) p26.push_back(p);
                if(p.pT >= 1.4 && p.pT < 1.6) p27.push_back(p);
                if(p.pT >= 1.6 && p.pT < 1.8) p28.push_back(p);
                if(p.pT >= 1.8 && p.pT < 2.0) p29.push_back(p);
            }

        }

        processEvent();

        p1.clear();
        p21.clear();
        p22.clear();
        p23.clear();
        p24.clear();
        p25.clear();
        p26.clear();
        p27.clear();
        p28.clear();
        p29.clear();

        if (!dataFile) break;
    }
}

void s1_parse_2partcorr()
{
    //Define 10 histograms
    for(int i=0; i<10; i++)
    {
        dhis.push_back(new TH1F(Form("d_%i",i),  "dhis", 50, -0.5*TMath::Pi(), 1.5*TMath::Pi()));
    }

    //Define count histgram
    hcount = new TH1F("hcount", "count", 6, 0, 6);

    //Make a file to store outputs
    TFile *fout = new TFile("out_two_particle.root","RECREATE");

    parseampt();

    //Save graphs
    for(int i=0; i<10; i++)
    {   
        dhis[i]->Write();
    }
    
    hcount->Write();
    fout->Close();
}


























