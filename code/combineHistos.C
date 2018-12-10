// Purpose: Combine different triggers into a single spectrum
// Author:  mikko.voutilainen@cern.ch
// Created: March 22, 2010
// Updated: June 2, 2015
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TObject.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile3D.h"

#include <iostream>
#include <map>
#include <string>

#include "settings.h"

using namespace std;

TH1D *recurseFile(TDirectory *indir, TDirectory *outdir, string hname = "hpt",
                  bool atBottom = false, TH1D *_hpt = 0, double etamid = 0);

// Map of pT ranges
std::map<std::string, std::pair<double, double>> _ptranges;

// Global variables (not pretty, but works)
TDirectory *_top = 0;

void combineHistos(std::string type) {

    bool _dt = (type == "DATA");
    bool _mc = !_dt;

    TDirectory *curdir = gDirectory;

    // Input file: normalized pT spectra
    TFile *fin = 
        new TFile(Form("../outputs/output-%s-2a.root", type.c_str()), "READ");
    assert(fin && !fin->IsZombie());
    
    // Top-level directory
    _top = gDirectory;

    // Output file: combined spectra 
    TFile *fout =
        new TFile(Form("../outputs/output-%s-2b.root", type.c_str()), "RECREATE");
    assert(fout && !fout->IsZombie());


    std::cout << "Calling combineHistos(" << type << ");" << std::endl;
    std::cout << "Input file " << fin->GetName() << std::endl;
    std::cout << "Output file " << fout->GetName() << std::endl;
    std::cout << "Starting recursions. These may take a few seconds" << std::endl;

    // Store pT ranges to a nice map
    for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {

        std::string name = _jp_triggers[itrg];
        double lower = _jp_trigranges[itrg][0];
        double upper = _jp_trigranges[itrg][1];

        _ptranges[name] = pair<double, double>(lower, upper);
    }


    // Loop over all the directories recursively
    recurseFile(fin, fout, "hpt");
    recurseFile(fin, fout, "hpt_pre");
    if (_dt)
        recurseFile(fin, fout, "hlumi");
    if (_mc)
        recurseFile(fin, fout, "hpt_g0tw");

    curdir->cd();

    std::cout << std::endl << "Output stored in " << fout->GetName() << std::endl;
    
    // Close files
    fout->Close();
    fout->Delete();
    std::cout << "Output file closed" << std::endl;

    fin->Close();
    fin->Delete();
    std::cout << "Input file closed" << std::endl;
}


TH1D *recurseFile(TDirectory *indir, TDirectory *outdir, string hname,
                  bool atBottom, TH1D *_hpt, double etamid) {
    TDirectory *curdir = gDirectory;

    // Recurse over all directories
    TList *keys = indir->GetListOfKeys();
    TListIter itkey(keys);
    TObject *key, *obj;

    while ((key = itkey.Next())) {

        std::string classname = ((TKey *)key)->GetClassName();
        std::string kname = key->GetName();

        // Search for a directory or histogram
        if (classname == "TDirectoryFile" || kname == hname) {
            obj = ((TKey *)key)->ReadObj();
            assert(obj);
        } 
        else {
            continue;
        }

        // Found a subdirectory: copy it to output and go deeper
        if (classname == "TDirectoryFile" && obj->InheritsFrom("TDirectory")) {
            TDirectory *outdir2 = outdir;

            if (!atBottom) {
                if (outdir->FindKey(obj->GetName()) == 0)
                    outdir->mkdir(obj->GetName());
            
                assert(outdir->cd(key->GetName()));
                outdir2 = outdir->GetDirectory(key->GetName());
            
                assert(outdir2);
                outdir2->cd();
            
                if (_debug)
                    cout << key->GetName() << endl;
            
            } else if (_debug)
                std::cout << key->GetName() << " (at bottom)" << endl;

            assert(indir->cd(obj->GetName()));
            
            TDirectory *indir2 = indir->GetDirectory(obj->GetName());
            assert(indir2);
            indir2->cd();

            // Check if directory name contains information on eta bin width
            // If yes, the next level is the bottom level with triggers
            // set flag and reset the combined histogram pointer
            float etamin, etamax;
            int valuesRead = sscanf(indir2->GetName(), "Eta_%f-%f", &etamin, &etamax);
            
            if (valuesRead == 2 && etamax > etamin) {

                float etamid = 0.5 * (etamin + etamax);
                _hpt = recurseFile(indir2, outdir2, hname, true, _hpt, etamid);

                if (_hpt) {
                    outdir2->cd();
                    _hpt->Write();
                    _hpt = 0;
                }
            } else {
                _hpt = recurseFile(indir2, outdir2, hname, atBottom, _hpt, etamid);
            }

        } 

        // Serch for a histogram
        else if (kname == hname && obj->InheritsFrom("TH1")) {
            
            TH1D *hpt = (TH1D *)obj;

            if (_hpt == 0) {
                outdir->cd();
                _hpt = (TH1D *)hpt->Clone(hpt->GetName());
                _hpt->Reset();
                indir->cd();
                assert(_hpt);

                if (_debug)
                    cout << "Cloned _" << hpt->GetName() << endl;
            }

            assert(_ptranges.find(indir->GetName()) != _ptranges.end() && "pT range not found for directory ");
            
            double ptmin = _ptranges[indir->GetName()].first;
            double ptmax = _ptranges[indir->GetName()].second;

            // Copy histogram values to output
            for (int i = 1; i != _hpt->GetNbinsX() + 1; ++i) {   
                double pt = _hpt->GetBinCenter(i);
                if (ptmin < pt && pt < ptmax) {
                    _hpt->SetBinContent(i, hpt->GetBinContent(i));
                    _hpt->SetBinError(i, hpt->GetBinError(i));
                } 
            }     
        } 
        // Free memory
        obj->Delete();
    } 

    // Free memory
    if (indir == _top) {
        std::cout << ".";;
    }

    curdir->cd();

    return _hpt;
} 
