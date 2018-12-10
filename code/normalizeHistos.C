// Purpose: Normalize inclusive jet analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: March 21, 2010
// Updated: June 2, 2015
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"

#include <iostream>

using namespace std;

#include "settings.h"

void recurseFile(TDirectory *indir, TDirectory *outdir,
         double etawid = 1., double etamid = 0.);

// Use this to fix luminosity
bool _mc = false;
bool _dt = true;

void normalizeHistos(string type) {

    assert(type==_jp_type);
    _mc = (type=="MC" || type=="HW");
    _dt = (type=="DATA");

    // Input file
    TFile *fin = new TFile(Form("../outputs/output-%s-1.root",type.c_str()),"READ");
    assert(fin && !fin->IsZombie());

    // Output file
    TFile *fout = new TFile(Form("../outputs/output-%s-2a.root",type.c_str()),"RECREATE");
    assert(fout && !fout->IsZombie());

    std::cout << "Calling normalizeHistos("<< type <<");" << std::endl;
    std::cout << "Input file " << fin->GetName() << std::endl;
    std::cout << "Output file " << fout->GetName() << std::endl;
    std::cout << "Starting recursive loop. This may take a minute" << std::endl << flush;

    // Loop over all the directories recursively
    recurseFile(fin, fout);

    std::cout << std::endl;
    std::cout << "Recursive loop done." << std::endl;
    std::cout << "Writing output to " << fout->GetName() << std::endl;
    std::cout << "This may again take a minute" << std::endl << flush;
    fout->Write();
    std::cout << "Output written in " << fout->GetName() << std::endl;
    fout->Close();
    std::cout << "Output file closed" << std::endl;
    fout->Delete();
    std::cout << "Output file pointer deleted" << std::endl << flush;


    fin->Close();
    fin->Delete();

} // normalizeHistos


void recurseFile(TDirectory *indir, TDirectory *outdir,
         double etawid, double etamid) {

    TDirectory *curdir = gDirectory;

    // Automatically go through the list of keys (directories)
    TList *keys = indir->GetListOfKeys();
    TListIter itkey(keys);
    TObject *key, *obj;
    TDirectory *dir;

    while ( (key = itkey.Next()) ) {

        if (_debug) 
        	std::cout << key->GetName() << std::endl << flush;
    
        obj = ((TKey*)key)->ReadObj(); assert(obj);
        dir = indir;
   
	string *fileName = new string(key->GetName());

        // Found a subdirectory: copy it to output and go deeper
        if (!(fileName->compare("MetSumetRatio") == 0) && obj->InheritsFrom("TDirectory")) {


            outdir->mkdir(obj->GetName());            
            assert(outdir->cd(obj->GetName()));

            TDirectory *outdir2 = outdir->GetDirectory(obj->GetName()); 
            assert(outdir2);
            
            outdir2->cd();
      
            assert(indir->cd(obj->GetName()));
            TDirectory *indir2 = indir->GetDirectory(obj->GetName());
            indir2->cd();

            // Check if directory name contains information on eta bin width
            float etamin, etamax;
            int valuesRead = sscanf(indir->GetName(),"Eta_%f-%f",&etamin,&etamax);
            if ( valuesRead == 2 && (etamax>etamin) ) {
                etawid = 2.*(etamax-etamin);
                etamid = 0.5*(etamax+etamin);
            }
    
            recurseFile(indir2, outdir2, etawid, etamid);
        } 

        // Found a plot: normalize if hpt, then copy to output
        else if (obj->InheritsFrom("TH1")) {

            outdir->cd();
            TObject *obj2 = obj->Clone(obj->GetName());

            // Normalize hpt 
            if (string(obj2->GetName())=="hpt"      || 
                string(obj2->GetName())=="hpt_pre"  ||
                string(obj2->GetName())=="hpt_g0tw"  ) {

                // Progress bar
                std::cout << ".";
            
                TH1D *hpt = (TH1D*)obj2;

                bool ispre = (TString(obj2->GetName()).Contains("_pre"));
            
                TH1D *hlumi = (TH1D*)dir->Get("hlumi"); 
                assert(hlumi && "Luminosity not found!");

                for (int i = 1; i != hpt->GetNbinsX()+1; ++i) {

                    // Normalization for bin width in y, pT
                    double norm = hpt->GetBinWidth(i) * etawid;
                    double pt = hpt->GetBinCenter(i);
                
                    // Luminosity normalization, except for prescaled histograms
                    if (hlumi->GetBinContent(i) != 0 && !ispre && _dt ) {
                        norm *= hlumi->GetBinContent(i);
                    }

                    assert(norm != 0 && "Invalid normalization!");

                    // Set values
                    hpt->SetBinContent(i, hpt->GetBinContent(i) / norm);
                    hpt->SetBinError(i, hpt->GetBinError(i) / norm);
            
                } 
            } 
            dir->cd();
        } 
    } 
    curdir->cd();
}
