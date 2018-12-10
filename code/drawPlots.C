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

#include "settings.h"
#include "tdrstyle_mod15.C"

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

        // Vector of eta ranges
	vector<string> etas;
 	etas.push_back("Eta_0.0-0.5");
  	etas.push_back("Eta_0.5-1.0");
  	etas.push_back("Eta_1.0-1.5");
  	etas.push_back("Eta_1.5-2.0");
  	etas.push_back("Eta_2.0-2.5");
  	etas.push_back("Eta_2.5-3.0");

        // Eta labels for the graph 
	map<string, string> label;
	label["Eta_0.0-0.5"] = "       |y|<0.5 (#times10^{5})";
 	label["Eta_0.5-1.0"] = "0.5#leq|y|<1.0 (#times10^{4})";
  	label["Eta_1.0-1.5"] = "1.0#leq|y|<1.5 (#times10^{3})";
  	label["Eta_1.5-2.0"] = "1.5#leq|y|<2.0 (#times10^{2})";
  	label["Eta_2.0-2.5"] = "2.0#leq|y|<2.5 (#times10)";
 	label["Eta_2.5-3.0"] = "2.5#leq|y|<3.0";

        // Markers for the graph
	map<string, int> marker;
	marker["Eta_0.0-0.5"] = kFullCircle;
	marker["Eta_0.5-1.0"] = kOpenCircle;
	marker["Eta_1.0-1.5"] = kFullSquare;
	marker["Eta_1.5-2.0"] = kOpenSquare;
	marker["Eta_2.0-2.5"] = kFullTriangleDown;
	marker["Eta_2.5-3.0"] = kOpenTriangleDown;

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

        // Vector of trigger names
	vector<string> trigger;
	trigger.push_back("jt370");
	trigger.push_back("jt240");
	trigger.push_back("jt190");
	trigger.push_back("jt110");
	trigger.push_back("jt60");
	trigger.push_back("jt30");

        // Trigger labels for the graph
	map<string, string> label;
	label["jt370"] = "Jet370";
	label["jt240"] = "Jet240";
	label["jt190"] = "Jet190";
	label["jt110"] = "Jet110";
	label["jt60"] = "Jet60";
	label["jt30"] = "Jet30";

        // Colour of the histograms
	map<string, int> colour;
	colour["jt370"] = kAzure-1;
	colour["jt240"] = kOrange+8;	
	colour["jt190"] = kViolet-4;
	colour["jt110"] = kTeal+7;
	colour["jt60"] = kPink-4;
	colour["jt30"] = kCyan+3;

        // Map of pT ranges for each trigger
        std::map<std::string, std::pair<double, double>> _ptranges;

        // Store pT ranges and their corresponding trigger in a map
        for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {

          std::string name = _jp_triggers[itrg];
          double lower = _jp_trigranges[itrg][0];
          double upper = _jp_trigranges[itrg][1];

          _ptranges[name] = pair<double, double>(lower, upper);
         }

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
		npt->SetMarkerStyle(kFullTriangleDown);
		npt->SetMarkerColor(colour[trigger[i]]);
		npt->SetMarkerSize(1.25);

		TH1D *on = (TH1D*)npt->Clone("on");

                // Find pT ranges for this trigger	
	        double ptmin = _ptranges[trigger[i]].first;
                double ptmax = _ptranges[trigger[i]].second;

                // Find bin corresponding to the pT value
                Int_t binmin = on->GetXaxis()->FindBin(ptmin);	
                Int_t binmax = on->GetXaxis()->FindBin(ptmax);	

		on->GetXaxis()->SetRange(binmin,binmax);

		leg->AddEntry(npt, label[trigger[i]].c_str(), "p");

		if (!strcmp(trigger[i].c_str(),"jt240") || !strcmp(trigger[i].c_str(),"jt190") || !strcmp(trigger[i].c_str(),"jt60")) {
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
