// Purpose:  Create basic histograms for inclusive jets analysis
//
// Author:   mikko.voutilainen@cern.ch
// Created:  March 20, 2010
// Updated:  June 8, 2015
#include "basicHistos.h"
#include "settings.h"

#include "TMath.h"

#include <iostream>

using namespace std;

basicHistos::basicHistos(TDirectory *dir, std::string trigname,
                         double ymin, double ymax,
                         double pttrg, double ptmin, double ptmax, bool ismc)
                        : lumsum(0) {
    
    TDirectory *curdir = gDirectory;
    dir->cd();
    this->dir = dir;

    // phase space
    this->trigname = trigname;
    this->ymin = ymin;
    this->ymax = ymax;
    this->pttrg = pttrg;
    this->ptmin = ptmin;
    this->ptmax = ptmax;
    this->ismc = ismc;

    // Once and for all (even if few too many with Sumw2)
    TH1::SetDefaultSumw2(kTRUE);

    // Binning agreed within JTF: pT>100 GeV from CaloJet resolutions,
    // pT<100 GeV to optimize bin widths for PFJets and b-tagging
    // (little higher than resolution, but fairly flat relative width)
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/QCDAnalysis/HighPtJetAnalysis/interface/DefaultPtBins.h?revision=1.2&view=markup
    const double x0[] = {
        1,    5,    6,    8,    10,   12,   15,   18,   21,   24,   28,   32,
        37,   43,   49,   56,   64,   74,   84,   97,   114,  133,  153,  174,
        196,  220,  245,  272,  300,  330,  362,  395,  430,  468,  507,  548,
        592,  638,  686,  737,  790,  846,  905,  967,  1032, 1101, 1172, 1248,
        1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
        2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
        4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};
    const int nx0 = sizeof(x0) / sizeof(x0[0]) - 1;
    //
    const double *bx0 = &x0[0];
    const int nbx0 = nx0;

    // Wider version of the binning for less statistical scatter for b-jets
    const double xW[] = {1,    5,    15,   24,   37,   56,   84,   114,
                         153,  196,  245,  330,  430,  548,  686,  846,
                         1032, 1248, 1497, 1784, 2116, 2500, 2941, 3450,
                         3637, 4252, 4961, 5777, 6717, 7000};
    const int nxW = sizeof(xW) / sizeof(xW[0]) - 1;
    //
    const double *bxW = &xW[0];
    const int nbxW = nxW;

    // Optimized binning created by optimizeBins.C ("MC"; lumi 1000/pb, eff
    // 1e+10%)
    // Using NLOxNP theory fit as input when available
    const int neta = 9;
    const int nbins = 65;
    double vx[neta][nbins] = {
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
         1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
         2640, 2787, 2941, 3103, 3273, 3450, 3832, 6076, 6389}, // Eta_0.0-0.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
         1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
         2640, 2787, 2941, 3103, 3273, 3637, 5220, 5492, 0}, // Eta_0.5-1.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
         1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
         2640, 2941, 3832, 4037, 0, 0, 0, 0, 0}, // Eta_1.0-1.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
         1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2500, 2640, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0}, // Eta_1.5-2.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
         1327, 1410, 1497, 1588, 1684, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0}, // Eta_2.0-2.5
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Eta_2.5-3.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, // Eta_2.5-3.0
        {10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114,
         133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507,
         548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}; // Eta_2.5-3.0

    // Rapidity bin index
    int ieta = int(0.5 * (ymin + ymax) / 0.5);
    assert(ieta < neta);

    // pT bins
    vector<double> x;
    for (int i = 0; i != nbins && vx[ieta][i] != 0; ++i) {
        x.push_back(vx[ieta][i]);
    } 

    // Number of pT bins
    const int nx = x.size() - 1;


    // Raw spectrum
    hpt = new TH1D("hpt", "", nx, &x[0]);
    hpt_pre = new TH1D("hpt_pre", "", nx, &x[0]); // prescale weighed

    // Luminosity
    hlumi = new TH1D("hlumi", "", nx, &x[0]);

    // Generator spectrum
    hpt_g0tw = new TH1D("hpt_g0tw", "", nx, &x[0]); // _mc per trigger

    curdir->cd();
}

basicHistos::~basicHistos() {
    dir->cd();
    dir->Write();
    delete dir;
};
