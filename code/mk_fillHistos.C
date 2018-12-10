// Purpose: Fill jet physics analysis histograms
// Author:  mikko.voutilainen@cern.ch
// Created: March 20, 2010


{
    // Compile code
    gROOT->ProcessLine(".L tools.C+");
    gROOT->ProcessLine(".L basicHistos.C+");

    gROOT->ProcessLine(".L fillHistos.C+g"); // +g for assert to work


    // Load tuples from file
    TChain *chain = new TChain("ak5ak7/OpenDataTree");
    assert(chain);

    if (_jp_type == "DATA") {
        std::cout << "Load trees..." << std::endl;
        
        chain->AddFile("root://eospublic.cern.ch//eos/opendata/cms/Run2012A/Jet/jettuples/OpenDataTuple-Data-Jet-Run2011A.root");

        std::cout << "Got " << chain->GetEntries() << " entries" << std::endl;
    }
    else if (_jp_type == "MC") 
    {
        std::cout << "Load trees..." << std::endl;

        chain->AddFile("root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/QCD_Pt-15to1000_TuneZ2_7TeV_pythia6/jettuples/OpenDataTuple-MC-QCD_Pt-15to1000_TuneZ2_7TeV_pythia6.root"); 

        std::cout << "Got " << chain->GetEntries() << " entries" << std::endl;
    }
    else {
        std::cout << "Invalid '_jp_type' value, take a look at settings.h!" << std::endl;
    }
    // Awkward patch for ROOT6:
    fillHistos(chain);
}
