// Purpose: Draw plots of analysis stages
// Author:  adelina.eleonora.lintuluoto@cern.ch
// Created: August 27, 2018
// Updated: 

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TMultiGraph.h"
#include "TMath.h"

//#include "tdrstyle_mod15.C"

#include <vector>

using namespace std;

// possible commands
void drawEtaSpectra();
void drawJetsPerBin();
void drawMetSumetRatio();

void drawPlots() {

	drawEtaSpectra();
	drawJetsPerBin();
        drawMetSumetRatio();
}

void drawEtaSpectra() {

	vector<string> etas;
 	etas.push_back("Eta_0.0-0.5");
  	etas.push_back("Eta_0.5-1.0");
  	etas.push_back("Eta_1.0-1.5");
  	etas.push_back("Eta_1.5-2.0");
  	etas.push_back("Eta_2.0-2.5");
  	etas.push_back("Eta_2.5-3.0");

	map<string, string> label;
	label["Eta_0.0-0.5"] = "       |y|<0.5 (#times10^{5})";
 	label["Eta_0.5-1.0"] = "0.5#leq|y|<1.0 (#times10^{4})";
  	label["Eta_1.0-1.5"] = "1.0#leq|y|<1.5 (#times10^{3})";
  	label["Eta_1.5-2.0"] = "1.5#leq|y|<2.0 (#times10^{2})";
  	label["Eta_2.0-2.5"] = "2.0#leq|y|<2.5 (#times10)";
 	label["Eta_2.5-3.0"] = "2.5#leq|y|<3.0";

	map<string, int> marker;
	marker["Eta_0.0-0.5"] = kFullCircle;
	marker["Eta_0.5-1.0"] = kOpenCircle;
	marker["Eta_1.0-1.5"] = kFullSquare;
	marker["Eta_1.5-2.0"] = kOpenSquare;
	marker["Eta_2.0-2.5"] = kFullTriangleDown;
	marker["Eta_2.5-3.0"] = kOpenTriangleDown;

	map<string, int> colour;
	colour["jt370"] = kAzure-1;
	colour["jt240"] = kOrange+8;	
	colour["jt190"] = kCyan+1;
	colour["jt150"] = kViolet-4;
	colour["jt110"] = kPink-4;
	colour["jt80"] = kGreen+2;
	
	double c = 10;
  	double scale[] = {pow(c,5),pow(c,4),pow(c,3),pow(c,2),pow(c,1),pow(c,0.)};

	// input file
	TFile *f = new TFile("../outputs/output-DATA-3.root","READ");
	assert (f && !f->IsZombie());

  	TDirectory *curdir = gDirectory;
	curdir->cd();

	double xmaxeta = 2500;
	double xmineta = 20;
  	TH1D *h = new TH1D("h","",int(xmaxeta-xmineta),xmineta,xmaxeta);
  	h->SetMinimum(1e-5);
  	h->SetMaximum(1e14);
  	h->GetYaxis()->SetTitle("d^{2}#sigma/dp_{T}dy (pb/GeV)");
  	h->GetXaxis()->SetTitle("Jet p_{T} (GeV)");
  	h->GetXaxis()->SetMoreLogLabels();
  	h->GetXaxis()->SetNoExponent();

	TF1 *f1 = new TF1("f1","[0]*pow(x,[1])"
		    "*pow(1-(2*1.*x*cosh([3]))/7000.,[2])",
		    0,2000);
  	f1->SetParameters(1.2e14,-5,10,2.5);
  	f1->SetParLimits(0, 1e11, 1e17);
  	f1->SetParLimits(1, -3, -7);
  	f1->SetParLimits(2, 8, 12);
  	f1->SetParLimits(3, 0, 0.5);
	f1->SetLineColorAlpha(kRed, 0.45);
	f1->SetLineStyle(1);

	TCanvas *c1 = tdrCanvas("c1",h,14,11,kSquare);

	TLegend *leg = tdrLeg(0.5886288,0.6515679,0.8896321,0.8815331);
  	leg->SetTextSize(0.035);

	TLegend *leg2 = tdrLeg(0.2123746,0.1114983,0.4130435,0.3118467);
	leg2->SetTextSize(0.035);

	TLatex *tex = new TLatex();
	tex->SetTextSize(0.035);
	tex->SetNDC();
	tex->DrawLatex(0.40,0.85,"Anti-k_{T} R=0.5");
	
	// values of eta
	double eta[] = {0, 0.5, 1, 1.5, 2, 2.5, 3};
	// parameter index zero, parameters corresponds to eta values {0.0-0.5, 0.5-1.0,.., 2.5-3.0}
	double par_0[] = {1e19, 1e18, 1e18, 1e16, 4e14, 5e13};
	// parameter limits, first pair belong to eta values 0.0-0.5, next pair for 0.5-1.0 .. etc. One pair = lower and upper limits
	double parlim_0[] = {1e16, 1e19, 1e15, 1e18, 1e15, 1e18, 1e15, 1e16, 1e13, 1e15, 1e13, 1e14}; 

	for (int i = 0; i != 6; i++) {
	
		// double differential sigma
		TGraph *g = (TGraph*)f->Get(Form("Standard/%s/gcorrpt",etas[i].c_str())); assert(g);

		for (int j = 0; j != g->GetN(); ++j) {
			// scale it according to 10^x[i] (x = {0,1,2,3,4,5})
			g->SetPoint(j, g->GetX()[j], g->GetY()[j] * scale[i]);		
		}

		g->SetMarkerStyle(marker[etas[i]]);
		g->SetMarkerSize(1.2);
		leg->AddEntry(g, label[etas[i]].c_str(), "p");
		g->Draw("x samep");
	
		// graph to fit to the data	
		TGraph *fit= (TGraph*)g->Clone("fit"); assert(fit);	
		
		f1->SetLineWidth(5);
		
		f1->FixParameter(0, par_0[i]);
		f1->SetParLimits(0, parlim_0[i], parlim_0[i+1]);	
		f1->FixParameter(3, (eta[i] - eta[i+1]) / 2);
		f1->SetParLimits(3, eta[i], eta[i+1]);
			
		if (i == 4) {

			f1->FixParameter(2, 4);
			f1->SetParLimits(2, 3, 4);				
			f1->FixParameter(3, 2.3);
			f1->SetParLimits(3, 2.2, 2.4);
		}

		if (i == 5) {

			f1->FixParameter(2,  4);
			f1->SetParLimits(2, 3, 4);				

		}

		int first = fit->GetX()[1];
		int last = fit->GetX()[fit->GetN()-1];
		
		if (i == 0) {
			fit->Fit("f1","","",first,1700);
		} else if (i < 3) {
			fit->Fit("f1","","",first,1500);
		} else if (i < 4 ) {
			fit->Fit("f1","","",first,1000);
		} else if (i < 5) {
			fit->Fit("f1","","",first,600);
		} else {
			fit->Fit("f1","","",first,400);
		}

		fit->FindObject("f1")->Draw("SAME");
		if (i == 0) leg2->AddEntry(f1, "Fit");
	} 

	gPad->SetLogy();
	gPad->SetLogx();
	gPad->RedrawAxis();	

	c1->SaveAs("../plots/eta_spectra.pdf");

} //drawEtaSpectra

void drawJetsPerBin() {

	bool _debug = false;

	vector<string> trigger;
	trigger.push_back("jt370");
	trigger.push_back("jt240");
	trigger.push_back("jt190");
	trigger.push_back("jt150");
	trigger.push_back("jt110");
	trigger.push_back("jt80");

	map<string, string> label;
	label["jt370"] = "Jet370";
	label["jt240"] = "Jet240";
	label["jt190"] = "Jet190";
	label["jt150"] = "Jet150";
	label["jt110"] = "Jet110";
	label["jt80"] = "Jet80";

	map<string, int> marker;
	marker["jt370"] = kOpenTriangleDown;
	marker["jt240"] = kFullTriangleDown;
	marker["jt190"] = kFullTriangleDown;
	marker["jt150"] = kOpenTriangleDown;
	marker["jt110"] = kFullTriangleDown;
	marker["jt80"] = kOpenTriangleDown;

	map<string, int> colour;
	colour["jt370"] = kAzure-1;
	colour["jt240"] = kOrange+8;	
	colour["jt190"] = kViolet-4;
	colour["jt150"] = kTeal+7;
	colour["jt110"] = kPink-4;
	colour["jt80"] = kCyan+3;

	// contains bins with largest y-value per trigger
	vector<int> range;

	// input file
	TFile *f = new TFile("../outputs/output-DATA-2a.root","READ");
	assert (f && !f->IsZombie());

  	TDirectory *curdir = gDirectory;
	curdir->cd();

	TH1D *h = new TH1D("h", ";Jet p_{T} (GeV);Number of jets per bin",100,0,2000);
	h->SetMaximum(450e3);
	h->SetMinimum(0);
	h->GetXaxis()->SetMoreLogLabels();
	h->GetXaxis()->SetNoExponent();

	TCanvas *c2 = tdrCanvas("c2",h,14,11,kSquare);

	TLegend *leg = tdrLeg(0.1839465,0.5400697,0.4849498,0.8205575);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);

	TLatex *tex = new TLatex();
  	tex->SetTextSize(0.04);
 	tex->SetNDC();
  	tex->DrawLatex(0.78,0.85,"#lbary#lbar < 0.5");

	for (int i = 0; i != 6; i++) {

		TH1D *pt = (TH1D*)f->Get(Form("Standard/Eta_0.0-0.5/%s/hpt",trigger[i].c_str())); assert(pt);

		TH1D *lumi = (TH1D*)f->Get(Form("Standard/Eta_0.0-0.5/%s/hlumi",trigger[i].c_str())); assert(lumi);

		TH1F *pre = (TH1F*)f->Get(Form("Standard/Eta_0.0-0.5/%s/hpt_pre",trigger[i].c_str())); assert(pre);
	
		TH1D *npt = (TH1D*)pt->Clone("npt"); assert(npt);

		TH1D *pt2 = (TH1D*)pt->Clone("pt2"); assert(pt2);

		for (int j = 1; j != npt->GetNbinsX() + 1; ++j) {

			npt->SetBinContent(j,npt->GetBinWidth(j)*npt->GetBinContent(j));
						
			npt->SetBinError(j,npt->GetBinWidth(j)*npt->GetBinError(j));

		}

		npt->Multiply(lumi);
		npt->SetMarkerStyle(marker[trigger[i]]);
		npt->SetMarkerColor(colour[trigger[i]]);
		npt->SetMarkerSize(1.25);

		TH1D *on = (TH1D*)npt->Clone("on");
		
		// number of bins
		int binsize = on->GetSize()-2;

		// add last bin
		if (i == 0) range.push_back(binsize);

		// find bin with highest y-value (omit bins before 13)
		int binmax = on->GetMaximumBin();
		range.insert(range.begin(),binmax);

		if (_debug) {

			cout << trigger[i] << endl;

			cout << Form("	first bin: %1.2d = %1.2f GeV, last bin: %1.2d = %1.2f GeV", range[0], on->GetXaxis()->GetBinCenter(range[0]), range[1], on->GetXaxis()->GetBinCenter(range[1])) << endl;
		
			cout << Form("	binsize: %1.2d, binmax: %1.2d, xvalue: %1.2f GeV, yvalue: %1.2f", binsize, binmax, on->GetXaxis()->GetBinCenter(binmax), on->GetBinContent(binmax)) << endl; 

		}
	
		on->GetXaxis()->SetRange(range[0],range[1]);

		leg->AddEntry(npt, label[trigger[i]].c_str(), "p");

		if (!strcmp(trigger[i].c_str(),"jt240") || !strcmp(trigger[i].c_str(),"jt190") || !strcmp(trigger[i].c_str(),"jt110")) {
			npt->SetFillColorAlpha(colour[trigger[i]], 0.35);
			npt->SetLineColorAlpha(colour[trigger[i]], 0.35);
			npt->Draw("hist same c");

			on->SetFillColorAlpha(colour[trigger[i]], 0.55);
			on->SetLineColorAlpha(colour[trigger[i]], 0.55);
			on->Draw("hist same c");
		} else {
			npt->SetLineColorAlpha(colour[trigger[i]], 0.35);
			npt->Draw("hist same cp");
			
		}
	
	} 

	gPad->SetLogx();
	gPad->RedrawAxis();

	c2->SaveAs("../plots/jets_per_bin.pdf");
		
}//drawJetsPerBin

void drawMetSumetRatio() {

        // input file
        TFile *f = new TFile("../outputs/output-DATA-1.root","READ");
        assert (f && !f->IsZombie());

        setTDRStyle();
        TDirectory *curdir = gDirectory;
        curdir->cd();

        TH1D *h = new TH1D("h","",50,0,1);
        h->GetXaxis()->SetTitle("met/sumet");
        h->SetMaximum(1e4);
        h->SetMinimum(0);

        TCanvas *c = tdrCanvas("c",h,14,11,kSquare);
        h->GetYaxis()->SetTitleOffset(1.6);

        TLegend *leg = tdrLeg(0.3386288,0.7090592,0.7884615,0.8815331);
        leg->SetTextSize(0.04);
        leg->SetTextFont(42);

        TH1D *leq100 = (TH1D*)f->Get("NoEventSelection/MetSumetRatio/leq100"); assert(leq100);

        TH1D *geq500 = (TH1D*)f->Get("NoEventSelection/MetSumetRatio/geq500"); assert(geq500);

        leq100->SetLineColor(kRed);
        geq500->SetLineColor(kBlue);

        geq500->Draw("hist same c");
        leq100->Draw("hist same c");
        leg->AddEntry(geq500,"jet_pt[0] greater than 500 GeV" , "l");
        leg->AddEntry(leq100, "jet_pt[0] less than 100 GeV", "l");

        gPad->RedrawAxis();

        c->SaveAs("../plots/met_sumet_ratio_of_different_pt_events.pdf");
} //drawMetSumetRatio
