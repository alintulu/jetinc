// Purpose: Reformat theory curves
// Author:  mikko.voutilainen@cern.ch
// Created: June 1, 2010
// Updated: June 2, 2015
//#include "ptresolution.h"
#include "tools.h"
#include "settings.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"

#include <string>
#include <fstream>
#include <sstream>

using namespace std;

// Return largest of three numbers
inline double max3(double a, double b, double c) {
    return max(max(a, b), c);
}
// Return largest of three numbers
inline double min3(double a, double b, double c) {
    return min(min(a, b), c);
}



// Load various theory predictions (FastNLO, LO MC)
void theoryBin(TDirectory *din, TDirectory *dth, TDirectory *dout);

// Also load old experimental data
void dataBin(TDirectory *din, TDirectory *dout);

void theory2(string type, string dir, TDirectory *fin, TDirectory *fmc, TDirectory *fout);

void theory(string type) {

    TDirectory *curdir = gDirectory;

    // TFile *fin = new TFile(Form("outputs/output-%s-3a.root",type.c_str()),"READ");
    TFile *fin = new TFile(Form("../outputs/output-%s-2b.root", type.c_str()), "READ");
    assert(fin && !fin->IsZombie());

    // LO Pythia MC prediction
    TFile *fmc = new TFile("../outputs/output-MC-2b.root", "READ");
    assert(fmc && !fmc->IsZombie());

    // Output
    TFile *fout = new TFile(Form("../outputs/output-%s-2c.root", type.c_str()), "RECREATE");
    assert(fout && !fout->IsZombie());

    theory2(type, "Standard", fin, fmc, fout);   
    curdir->cd();
    theory2(type, "NoEventSelection", fin, fmc, fout);   

    fout->Write();
    fout->Close();
    fout->Delete();

    fin->Close();
    fin->Delete();

    fmc->Close();
}

void theory2(string type, string dir, TDirectory *fin, TDirectory *fmc, TDirectory *fout) {

    // Select top category
    fin->cd(dir.c_str());
    TDirectory *din0 = gDirectory;

    fmc->cd(dir.c_str());
    TDirectory *dmc0 = gDirectory;

    fout->mkdir(dir.c_str());
    fout->cd(dir.c_str());
    TDirectory *dout0 = gDirectory;

    // Automatically go through the list of keys (directories)
    TList *keys = din0->GetListOfKeys();
    TListIter itkey(keys);
    TObject *key, *obj;

    while ((key = itkey.Next())) {
        obj = ((TKey *)key)->ReadObj();
        assert(obj);

        // Found a subdirectory
        if (obj->InheritsFrom("TDirectory") &&
            string(obj->GetName()) != "Eta_0.0-1.3" &&
            string(obj->GetName()) != "Eta_3.0-3.2" &&
            string(obj->GetName()) != "Eta_3.2-4.7") {

            din0->cd(obj->GetName());
            TDirectory *din = gDirectory;

            dmc0->cd(obj->GetName());
            TDirectory *dmc = gDirectory;

            dout0->mkdir(obj->GetName());
            dout0->cd(obj->GetName());
            TDirectory *dout = gDirectory;


            // Process rapidity bin
            theoryBin(din, dmc, dout);

        } // inherits TDirectory
    }     // while

    cout << "Output stored in " << fout->GetName() << endl;

} // theory


void theoryBin(TDirectory *din, TDirectory *dth, TDirectory *dout) {

    float etamin = 0;
    float etamax = 0;
    sscanf(din->GetName(), "Eta_%f-%f", &etamin, &etamax);

    std::cout << etamin << ' ' << etamax << std::endl;

    int ieta = int((0.5 * (etamin + etamax)) / 0.5);
    int jeta = (_jp_algo == "AK7" ? min(4, ieta) : ieta);

    // inclusive jets
    TH1D *hpt = (TH1D *)din->Get("hpt");

    // theory curves
    // http://www-ekp.physik.uni-karlsruhe.de/~rabbertz/fastNLO_LHC/InclusiveJets/
    // -> fnl2342_cteq66_aspdf_full.root
    // Numbering scheme explained in
    // https://twiki.cern.ch/twiki/bin/view/CMS/CMSfastNLO
    // 2-point to 6-point theory uncertainty:
    // h200X00->h300X09, h400X00->h300X08
    TH1D *hnlo(0); 

    cout << din->GetName() << " ieta=" << ieta << " jeta=" << jeta << endl;

    TH1D *hmc = (TH1D *)dth->Get("hpt_g0tw");
    assert(hmc);

    // Make sure new histograms get created in the output file
    dout->cd();

    // Initial fit of the NLO curve
    TF1 *fnlo = new TF1("fnlo", "[0]*exp([1]/x)*pow(x,[2])*pow(1-x*cosh([4])/3500.,[3])", 10., 1000.);
    fnlo->SetParameters(2e14, -18, -5.2, 8.9, 0.5 * ieta);
    fnlo->FixParameter(4, 0.5 * ieta);

    hnlo = hmc;
    hnlo->Fit(fnlo, "QRN");

    // Graph of theory points with centered bins
    const double minerr = 0.02;
    TGraphErrors *gnlo = new TGraphErrors(0);
    TGraphErrors *gnlo2 = new TGraphErrors(0);   // above + minerr
    TGraphErrors *gnlocut = new TGraphErrors(0); // only up to expected pT
    gnlo->SetName("gnlo");
    gnlo2->SetName("gnlo2");
    gnlocut->SetName("gnlocut");

    for (int i = 1; i != hnlo->GetNbinsX() + 1; ++i) {
        double y = hnlo->GetBinContent(i);
        double dy = hnlo->GetBinError(i);
        double ptmin = hnlo->GetBinLowEdge(i);
        double ptmax = hnlo->GetBinLowEdge(i + 1);

        double y0 = fnlo->Integral(ptmin, ptmax) / (ptmax - ptmin);
        double x = fnlo->GetX(y0, ptmin, ptmax);

        int n = gnlo->GetN();
        tools::SetPoint(gnlo, n, x, y, 0, dy);
        tools::SetPoint(gnlo2, n, x, y, 0, tools::oplus(dy, minerr * y));
        if (y * (ptmax - ptmin) > 0.0001) { // for 10/fb
            int m = gnlocut->GetN();
            tools::SetPoint(gnlocut, m, x, y, 0, 0);
        }
    }

    gnlo2->Fit(fnlo, "QRN");

    // Divide graph with fit to check stability
    TGraphErrors *gnlofit = new TGraphErrors(0);
    gnlofit->SetName("gnlofit");
    for (int i = 0; i != gnlo->GetN(); ++i) {
        double x, y, ex, ey;
        tools::GetPoint(gnlo, i, x, y, ex, ey);
        double f = fnlo->Eval(x);
        tools::SetPoint(gnlofit, i, x, y / f, ex, ey / f);
    }

    // Rebin theory to match data bins
    TH1D *hnlo0 = hnlo;
    hnlo0->SetName("hnlo0");
    hnlo = tools::Rebin(hnlo0, hpt);
    hnlo->SetName("hnlo");


    dout->cd();
    gnlo->Write();
    gnlocut->Write();
    gnlofit->Write();
    fnlo->Write();

    din->cd();
} // theoryBin
