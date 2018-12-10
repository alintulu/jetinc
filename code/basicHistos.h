// Purpose:  Define basic histograms for jet physics studies
// Author:   mikko.voutilainen@cern.ch
// Created:  March 20, 2010
// Updated:  June 9, 2015


#ifndef __basicHistos_h__
#define __basicHistos_h__

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TDirectory.h"

#include <string>
#include <map>
#include <vector>

class basicHistos {

 public:

  // Phase space
  std::string trigname;
  double ymin;
  double ymax;
  double pttrg;
  double ptmin;
  double ptmax;
  bool ismc;

  // Luminosity
  TH1D *hlumi;
  std::map<int, std::map<int, float> > lums;
  double lumsum;

  // Raw spectra
  TH1D *hpt;
  TH1D *hpt_pre;

  // Unbiased generator spectrum
  TH1D *hpt_g0tw;

  basicHistos(TDirectory *dir, std::string trigname="", 
	            double ymin = 0., double ymax = 2.0,
	            double pttrg = 10., double ptmin = 10., double ptmax = 50.,
	            bool ismc = false);
  ~basicHistos();

 private:
  TDirectory *dir;
};

#endif
