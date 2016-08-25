//////////////////////////////////////////////////////////////////////////////
// Authors:     Alfredo Gurrola, Andres Florez                              //
// contact:                                                                 //
//   Alfredo.Gurrola@cern.ch       (Vanderbilt University)                  //
//   Andres.Florez@cern.ch         (Los Andes University)                   //
//////////////////////////////////////////////////////////////////////////////

#ifndef BSM3GNormalizer_h
#define BSM3GNormalizer_h

// system include files
#include <memory>

// user include files
#include <Math/VectorUtil.h>
#include <fstream>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TFile.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <utility>
#include <TROOT.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TLorentzVector.h>
#include <TEnv.h>
#include <TError.h>
#include <TCollection.h>
#include <TKey.h>
#include <TGraphAsymmErrors.h>
#include <TClass.h>

using namespace std;

class BSM3GNormalizer {
public:
  BSM3GNormalizer(char*);
  ~BSM3GNormalizer();

private:

  virtual void beginJob();
  virtual void endJob();
  virtual void getInputs(char*);
  virtual void createLogFile();
  virtual void grabYieldsANDgrabHistos();
  virtual void calculateEfficienciesAndErrors();
  virtual void NormalizeHistos();
  virtual void CreateProbHistos();

  // initialize variables
  float effSkim;
  vector<string> inRootFiles;
  vector<string> inDirectories;
  vector<string> inProcess;
  vector<float> inScaleFactor;
  vector<float> inScaleFactorError;
  vector<float> x_section;
  vector<float> theCumulativeEfficiency;
  vector<float> RelativeEffVector;
  vector<float> RelativeEffErrorVector;
  vector<float> TotalEffVector;
  vector<float> TotalEffErrorVector;
  vector<float> theEventsAnalyzed;
  vector<float> theEventsPassing;
  vector<string> theHistNameStrings;
  vector<double> maxHistContents;
  vector<TH1F*> HistList;
  vector<TH1F*> HistList2;
//  TList *HistList = new TList();
//  TList *HistList2 = new TList();
  int nHistList;
  int nHistList2;
  int nHistos;
  string outputLogFile;
  string IsData;
  float MaxEventsAnalyzed;
  string outputRootFileForNormalizedHistos;
  string outputRootFileForProbabilityHistos;
  float lumi;
  string lumi_string;
  TFile *theCurrentFile;

};
#endif
