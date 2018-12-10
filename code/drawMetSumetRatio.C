#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "TMultiGraph.h"
#include "TMath.h"

#include "tools.C"
#include "tdrstyle_mod15.C"

#include <vector>

using namespace std;

void drawMetSumetRatio() {

	// input file
	TFile *f = new TFile("output-DATA-1.root","READ");
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

	c->SaveAs("met_sumet_ratio_of_different_pt_events.pdf");
}
