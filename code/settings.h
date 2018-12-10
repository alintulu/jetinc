// Purpose: Analysis settings for jet physics package
// Author:  mikko.voutilainen@cern.ch
// Created: June 1, 2015
// Modified: August 23, 2016 
#ifndef __settings_h__
#define __settings_h__
#include <string>
#include <iostream>

// All the settings variables are in global name space
// To avoid conflicts, and to make them easier to find
// global variables all start with _jp_ (for JetPhysics)

// Print out debugging information (for this only, no _jp_)
const bool _debug = false;

////////////////////////////////////////
// Necessary analysis settings/inputs //
////////////////////////////////////////

// Algorithm to use ("AK4", "AK5", "AK8")
std::string _jp_algo = "AK5";

// Data type ("DATA", "MC", or "HW")
std::string _jp_type = "MC";

// Number of events to process (-1 for all)
Long64_t _jp_nentries = -1;

// Calculate luminosity on the fly based on .csv file
const bool _jp_dolumi = true;
std::string _jp_lumifile = "lumicalc/pixellumi_by_LS.csv";

// List of triggers used in the analysis
const int _jp_ntrigger = 6;
std::string _jp_triggers[_jp_ntrigger] = {"jt30",  "jt60",    "jt110",
                                          "jt190", "jt240", "jt370"};

// Trigger pT thresholds and ranges
double _jp_trigthr[_jp_ntrigger]        = {       30,         60,         110,    
                                                 190,         240,        370};
double _jp_trigranges[_jp_ntrigger][2]  = {{0,   97},  {97, 174},   {174,  300},  
                                          {300, 362},  {362, 507}, {507, 3500}};

// Produce basic set of histograms
const bool _jp_doBasicHistos = true;

// Process pThatbins instead of flat sample
const bool _jp_pthatbins = true;
// For creating smearing matrix
const bool _jp_doMatrix = false;

// Only load selected branches 
// (significant speedup, but remember to enable all the right branches!)
const bool _jp_quick = true;

// Minimum and maximum pT range to be plotted and fitted
const double _jp_recopt = 24; 
const double _jp_fitptmin = 43;
const double _jp_xminpas = 56;
const double _jp_xmin = 24.;
const double _jp_xmax = 1999.;

// Center of mass energy
const double _jp_sqrts = 13000.; // GeV
const double _jp_emax = _jp_sqrts/2.; // GeV

#endif // __settings_h__
