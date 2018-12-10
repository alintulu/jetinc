#define fillHistos_cxx
#include "fillHistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void fillHistos::Loop()
{

    // Check that chain is valid
    if (fChain == 0) {
        return;
    }

    // Read recorded luminosities from text file
    if (_dt && _jp_dolumi) loadLumi(_jp_lumifile);

    // Number of events
    Long64_t nentries = 0; 
    if (_jp_nentries >= 0) {
        nentries = min(fChain->GetEntries(), _jp_nentries);
    }
    else {
        nentries = fChain->GetEntries();
    }
    
    // Switch all branches OFF
    fChain->SetBranchStatus("*", 0);

    fChain->SetBranchStatus("njet", 1);
    fChain->SetBranchStatus("jet_pt", 1);
    fChain->SetBranchStatus("jet_eta", 1);
    fChain->SetBranchStatus("jet_phi", 1);
    fChain->SetBranchStatus("jet_E", 1);

    fChain->SetBranchStatus("njet", 1);
    fChain->SetBranchStatus("jet_*", 1);

    // Switch on Monte Carlo branches
    if (_jp_type == "MC") {
        fChain->SetBranchStatus("ngen", 1);
        fChain->SetBranchStatus("gen_*", 1);
        fChain->SetBranchStatus("mcweight", 1);
    }

    fChain->SetBranchStatus("triggers", 1);
    fChain->SetBranchStatus("triggernames", 1);
    fChain->SetBranchStatus("prescales", 1);

    fChain->SetBranchStatus("run", 1);
    fChain->SetBranchStatus("lumi", 1);

    if (_jp_doBasicHistos) {
        initBasics("Standard");
    }
    
    // Event loop
    for (Long64_t jentry=0; jentry < nentries;jentry++) {

        // Error check
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        
        // Read event
        fChain->GetEntry(jentry);

        // Write this event to histograms
        if (_jp_doBasicHistos) {
            fillBasics("Standard");
        }    
    }
        
    // Delete histograms and close file properly
    if (_jp_doBasicHistos) {
        writeBasics(); 
    }
   
}


void fillHistos::initBasics(std::string name) {
// Initialize basic histograms for trigger and eta bins

    // Save current directory
    TDirectory *curdir = gDirectory;

    // Create output file
    TFile *f = (_outfile ? _outfile: new TFile(Form("output-%s-1.root", _type.c_str()), "RECREATE"));
    
    // Safety check
    assert(f && !f->IsZombie() && "Error while creating output file!");
    
    // Create top-level directory
    f->mkdir(name.c_str());
    assert(f->cd(name.c_str()) && "Error while creating directory");

    // Move into the directory
    TDirectory *topdir = f->GetDirectory(name.c_str());
    assert(topdir);    
    topdir->cd();

    // Rapidity bins + HF + barrel
    double y[] = {0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 4.7, 0., 1.3};
    const int ny = sizeof(y) / sizeof(y[0]) - 1;


    // Trigger bins
    // Efficient trigger pT ranges
    map< std::string, pair< double, double > > pt;
    for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {
        std::string trg = _jp_triggers[itrg];       // Trigger name
        double pt1      = _jp_trigranges[itrg][0];  // Lower bound
        double pt2      = _jp_trigranges[itrg][1];  // Upper bound

        pt[trg] = pair< double, double >(pt1, pt2);
    }

    // Trigger thresholds
    map< std::string, double > pttrg;
    for (int itrg = 0; itrg != _jp_ntrigger; ++itrg) {
        std::string trg = _jp_triggers[itrg]; // Trigger name
        double pt0      = _jp_trigthr[itrg];  // Threshold
        
        pttrg[trg] = pt0;
    }


    // Read first entry from the TTree to get trigger names 
    fChain->GetEntry(0);
    std::vector<std::string> trg_names = *triggernames;

    // Loop over rapidity and trigger bins
    for (int i = 0; i != ny; ++i) {
        if (y[i + 1] > y[i]) { // Create only real bins

            // Subdirectory for rapidity bin
            const char *yname = Form("Eta_%1.1f-%1.1f", y[i], y[i + 1]);
            assert(topdir);

            // Change to the subdirectory
            topdir->mkdir(yname);
            assert(topdir->cd(yname) && "Error while creating directory");
            
            TDirectory *ydir = topdir->GetDirectory(yname);
            assert(ydir);
            
            ydir->cd();
            
            // Loop over triggers
            for (std::string trg_name: trg_names) {

                // Trigger subdirectory
                const char *trg = trg_name.c_str();
                assert(ydir);
                
                // Create directory
                ydir->mkdir(trg);
                ydir->cd(trg);
                
                // Change directory
                TDirectory *dir = ydir->GetDirectory(trg);
                assert(dir);
                dir->cd();

                // Initialize histograms of this bin 
                basicHistos *h = new basicHistos(dir, trg, y[i], y[i + 1], pttrg[trg],
                                                 pt[trg].first, pt[trg].second, _mc);
            
                _histos[name].push_back(h);
            }
        }    
    }        

    _outfile = f;
    curdir->cd();

} 


// Loop over the histograms of the container and fill them
void fillHistos::fillBasics(std::string name) {
    
    for (auto const& h: _histos[name]) {
        fillBasic(h);
    }
}

// Fill basic histograms after applying pT and y cuts
void fillHistos::fillBasic(basicHistos *h) {

    // Check valididty
    assert(h && "Invalid set of histograms!");

    // Index of the trigger corresponding to this histogram
    unsigned int trg_index = find(triggernames->begin(), triggernames->end(), h->trigname) - triggernames->begin();
    
    // Return if not triggered
    bool fired = (triggers[trg_index]);
    if (!fired) {
        return;
    }



    // Check for missing prescale
    unsigned int prescale = prescales[trg_index];
    if (prescale == 0 && fired) {
        *ferr   << "Prescale zero for trigger " << h->trigname
                << " in run " << run << "!" << endl << flush; 
        assert(false);
    }


    // Luminosity information
    if (_dt && h->lums[run][lumi] == 0) {
        // Recorded luminosity in this lumi section
        double lum = _lums[run][lumi];
        
        // Count this lumisection in the sum (with prescale)
        h->lumsum += lum / prescale;

        // No double counting
        h->lums[run][lumi] = 1;
    }


    // Loop over jets of this event
    for (unsigned int i = 0; i != njet; ++i) {
        
        if (_debug) {
           std::cout << "Loop over jet " << i << "/" << njet << endl;
        }

        // Corrected 4-momentum
        double pt       = jet_pt[i];
        double eta      = jet_eta[i];
        double phi      = jet_phi[i];
        double energy   = jet_E[i];

        // Compute rapidity and save it
        p4.SetPtEtaPhiE(pt, eta, phi, energy);
        double y = p4.Rapidity();
        jet_y[i] = y; 


        // Absolute rapidity
        double y_abs = fabs(jet_y[i]);

        // Test whether jet belongs to this rapidity bin
        if (pt > _jp_recopt && h->ymin <= y_abs && y_abs < h->ymax)  {


            // Fill raw pT spectrum
            assert(h->hpt);
            if (_dt) {
                h->hpt->Fill(pt);
            }
            else if (_mc) {
                h->hpt->Fill(pt, mcweight);
            }

            // Fill prescaled pT spectrum
            if (_dt) {
                h->hpt_pre->Fill(pt, prescales[trg_index]);
            }

        }
    }  

    // Unbiased generator spectrum (needed for unfolding)
    if (_mc) {
        for (unsigned int i = 0; i != ngen; ++i) {
        
            // GenJet rapidity
            p4gen.SetPtEtaPhiE(gen_pt[i], gen_eta[i], gen_phi[i], gen_E[i]);
            double y = p4.Rapidity();

            if (h->ymin <= fabs(y) && fabs(y) < h->ymax) {
                h->hpt_g0tw->Fill(gen_pt[i], mcweight);

            }
        }
    }
    
} 

void fillHistos::writeBasics() {

    for (auto it = _histos.begin(); it != _histos.end(); ++it) {
        std::vector<basicHistos*> histoList = it->second;

        // Loop over bins of the output file
        for ( basicHistos *h : histoList) {

            // For data, save effective luminosity
            for (int j = 0; j != h->hlumi->GetNbinsX() + 1; ++j) {
                h->hlumi->SetBinContent(j, _dt ? h->lumsum : 1.);
            }

            // Delete histograms from memory
            delete h; 
        }             
    }   

    std::cout << "\nOutput stored in " << _outfile->GetName() << endl;

    _outfile->Close();
    _outfile->Delete();
}


// Load luminosity information
void fillHistos::loadLumi(const std::string filename) {

    // Open luminosity file
    std::ifstream CSVfile(filename);
    assert(CSVfile.is_open() && "Error while opening luminosity file!");

    // Read header
    std::string header;
    getline(CSVfile, header, '\r');
    std::cout << "CSV header: " << header << std::endl;
    
    // Require proper header format
    assert(header == string("Run:Fill,LS,UTCTime,Beam Status,E(GeV),Delivered(/ub),Recorded(/ub),avgPU"));
    
    // Total number of lumisections
    int nls = 0;

    // Total recorded luminosity
    double lumsum = 0;

    // One line of the file
    std::string line;
    while ( getline(CSVfile, line, '\r') ) {

        double rec;      // Recorded lumi in microbarns (ub)        
        unsigned int rn; // Run number
        unsigned int ls; // Lumisection number

        // Parse one line
        sscanf(line.c_str(),"%d:%*d,%d:%*d,%*d/%*d/%*d %*d:%*d:%*d,STABLE BEAMS,%*f,%*f,%lf,%*f", &rn, &ls, &rec);
        
        double lum = rec*1e-6;  // Convert ub to pb
        _lums[rn][ls] = lum;    // Save lumi section luminosity
        lumsum += lum;          // Add to luminosity sum 

        ++nls;
	assert(nls < 10e8 && "Error while reading luminosity info!");  
    }

    // Summary
    std::cout << "Called loadLumi(\"" << filename << "\"):" << endl;
    std::cout << "Loaded " << _lums.size() << " runs with "
              << nls << " lumi sections containing "
              << lumsum << " pb-1 of data " << std::endl;
    std::cout << "This corresponds to " << nls*23.3/3600
              << " hours of data-taking" << endl;
} 
