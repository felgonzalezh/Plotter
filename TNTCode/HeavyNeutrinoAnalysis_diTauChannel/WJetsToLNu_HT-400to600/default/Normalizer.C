#include <iostream>
#include <iomanip>
#include <vector>



void Normalizer() {

  #include <iostream>
  #include <iomanip>
  #include <vector>
  //gROOT->Reset();
  //  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  // open input file containing all information: root files, scale factors, cross-sections, etc
  ifstream inFile;
  char inputFilename[] = "Normalizer.in";
  inFile.open(inputFilename, ios::in);
  // if can't open input file, exit the code
  if (!inFile) {
    cerr << "Can't open input file " << inputFilename << endl;
    exit(1);
  }

  // initialize variables
  float effSkim;
  vector<string> inRootFiles;
  inRootFiles.clear();
  vector<string> inDirectories;
  inDirectories.clear();
  vector<string> inProcess;
  inProcess.clear();
  vector<float> inScaleFactor;
  inScaleFactor.clear();
  vector<float> inScaleFactorError;
  inScaleFactorError.clear();
  vector<float> x_section;
  x_section.clear();
  vector<float> theCumulativeEfficiency;
  theCumulativeEfficiency.clear();
  vector<float> RelativeEffVector;
  RelativeEffVector.clear();
  vector<float> RelativeEffErrorVector;
  RelativeEffErrorVector.clear();
  vector<float> theEventsAnalyzed;
  theEventsAnalyzed.clear();
  vector<float> theEventsPassing;
  theEventsPassing.clear();
  vector<string> theHistNameStrings;
  theHistNameStrings.clear();
  string inputString;
  string inputType;

  // grab all relevant information from the input file
  while (inFile >> inputType >> inputString) {
    if(inputType == "rootfile") {
      inRootFiles.push_back(inputString); // create a vector of root files containing histograms
    } else if(inputType == "directory") {
      inDirectories.push_back(inputString); // directory path within root file that contains the histograms
    } else if(inputType == "process") {
      inProcess.push_back(inputString); // name of the cut ... will be used to rename histograms
    } else if(inputType == "outputPostscriptFile") {
      string psFile1 = inputString; // name of output postscript file
    } else if(inputType == "outputRootFileForNormalizedHistos") {
      string outputRootFileForNormalizedHistos = inputString; // name of output root file containing normalized histograms
    } else if(inputType == "outputRootFileForProbabilityHistos") {
      string outputRootFileForProbabilityHistos = inputString; // name of output root file containing probability histograms
    } else if(inputType == "outputLogFile") {
      string outputLogFile = inputString; // name of output log file (e.g. contains efficiency table)
    } else if(inputType == "scaleFactor") {
      inScaleFactor.push_back(atof(inputString.c_str())); // scale factors used to recalculate efficiencies for a each cut
    } else if(inputType == "scaleFactorError") {
      inScaleFactorError.push_back(atof(inputString.c_str())); // scale factor error used to recalculate efficiency for each cut
    } else if(inputType == "luminosity") {
      float lumi = atof(inputString.c_str()); // luminosity
      string lumi_string = inputString;
    } else if(inputType == "effectiveXsection") {
      x_section.push_back(atof(inputString.c_str())); // cross-section
    } else if(inputType == "skimmingEff") {
      effSkim = atof(inputString.c_str()); // skimming efficiency
    } else if(inputType == "IsData") {
      string IsData = inputString;
    } else {
      cerr << "Incorrect input type " << inputType << endl; // exit code if unwanted input is specified
      exit(1);
    }
  }

  // create an output log file that will contain the cut flow efficiency table
  ofstream outFile;
  char outputFilename[] = outputLogFile.c_str();
  outFile.open(outputFilename, ios::out);
  // if output log file cannot be opened, exit the code
  if (!outFile) {
    cerr << "Can't open output file " << outputFilename << endl;
    exit(1);
  } else {
    outFile << "" << endl;
    outFile << "The following input was used: " << endl;
    outFile << "" << endl;
    outFile << "" << endl;
  }

  // strings needed to correctly write to a postscript (used later in the code)
  string psFile2 = psFile1+'('; // used if saving the first object to the postscript
  string psFile3 = psFile1+')'; // used if saving the last object to the postscript

  // initialization of variables
  vector<double> maxHistContents;
  maxHistContents.clear();
  HistList = new TList();
  int nHistList=0;
  HistList2 = new TList();
  int nHistList2=0;

  // loop over root files
  for(int j=0;j<inRootFiles.size();j++) {
    // printout information to log file
    outFile << "Name of Root File #" << (j+1) << " : " << inRootFiles.at(j) << endl;
    outFile << "Name of Directory #" << (j+1) << " : " << inDirectories.at(j) << endl;
    outFile << "Name of Process #" << (j+1) << "   : " << inProcess.at(j) << endl;
    outFile << "" << endl;
    TFile *example1 = new TFile(inRootFiles.at(j).c_str()); // open root file
    example1.cd(inDirectories.at(j).c_str()); // cd to appropriate directory
    TDirectory *current_sourcedir = gDirectory;
    TIter nextkey( current_sourcedir->GetListOfKeys() );
    TKey *key;
    int num=0;
    // loop over keys(ROOT objects) within the root file
    while((key = (TKey*)nextkey())) {
      TObject *obj = key->ReadObj();
      // only consider 1D histograms within the root file
      if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
        string histname = obj->GetName();
        TH1 *hobj = (TH1*)obj;
        // use the "Events" histogram to calculate cumulative efficiencies
        if( (hobj->GetYaxis()->GetNbins() == 1) && (histname == "Events_0") ) {
          // "Events" histogram contains 2 filled bins: bin 1 is the # of events analyzed; bin 2 is the # passing specified cuts
          if(hobj->GetBinContent(1) == 0) {
            cerr << "'Events' histogram contains zero entries in bin 1: 0 events were analyzed ..." << endl;
            exit(1);
          } else {
            /* calculate cumulative efficiency (at this stage in the code, the calculation is not complete because
               it does not incorporate skimming efficiencies */ 
            if(IsData == "1") {
              theEventsAnalyzed.push_back((double)hobj->GetBinContent(1)); // denominator for calculation of cumulative eff.
              theEventsPassing.push_back((double)hobj->GetBinContent(2)); // numerator for calculation of cumulative eff.
              theCumulativeEfficiency.push_back(((double)hobj->GetBinContent(2)) / ((double)hobj->GetBinContent(1))); // cumulative eff.
            } else {
//              theEventsAnalyzed.push_back((double)hobj->GetBinContent(3)); // denominator for calculation of cumulative eff.
              theEventsAnalyzed.push_back((double)hobj->GetBinContent(1)); // denominator for calculation of cumulative eff.
              theEventsPassing.push_back((double)hobj->GetBinContent(2)); // numerator for calculation of cumulative eff.
//              theCumulativeEfficiency.push_back(((double)hobj->GetBinContent(2)) / ((double)hobj->GetBinContent(3))); // cumulative eff.
              theCumulativeEfficiency.push_back(((double)hobj->GetBinContent(2)) / ((double)hobj->GetBinContent(1))); // cumulative eff.
//              effSkim = effSkim * (double)hobj->GetBinContent(3) / (double)hobj->GetBinContent(1);
            }
          }
        }
        // grab all 1D histograms within the root file and create a TList (i.e. vector) containing them to be used later
        if( (hobj->GetYaxis()->GetNbins() == 1) && (histname != "Events_0") ) {
          // rename histogram
          if(j < (inRootFiles.size() - 1)) {string histoName = histname + "_" + "After" + inProcess.at(j) + "Before" + inProcess.at(j+1);}
          else {string histoName = histname + "_" + "After" + inProcess.at(j);}
          string histoNameName = histoName + "2";
          h1 = new TH1F(histoName.c_str(),histoName.c_str(),hobj->GetXaxis()->GetNbins(),
                        hobj->GetXaxis()->GetXmin(),hobj->GetXaxis()->GetXmax());
          for(int i=0; i<=(h1->GetXaxis()->GetNbins() + 1); i++) {h1->SetBinContent(i,hobj->GetBinContent(i));}
          h2 = new TH1F(histoNameName.c_str(),histoNameName.c_str(),hobj->GetXaxis()->GetNbins(),
                        hobj->GetXaxis()->GetXmin(),hobj->GetXaxis()->GetXmax());
          for(int i=0; i<=(h2->GetXaxis()->GetNbins() + 1); i++) {h2->SetBinContent(i,hobj->GetBinContent(i));}
          HistList->Add(h1);
          nHistList++;
          HistList2->Add(h2);
          nHistList2++;
          theHistNameStrings.push_back(histname); // store default names - used later in the code
          double max = 0;
          for(int i=0; i<=(h1->GetXaxis()->GetNbins() + 1); i++) { if(h1->GetBinContent(i) > max) {max = h1->GetBinContent(i);} }
          if(j==0) {maxHistContents.push_back(max);}
          else {if(max > maxHistContents.at(num)) {maxHistContents.at(num) = max;}} // determine bin with largest entries
          num++;
        }
      }
    }
    // make sure that each root file contains the same histograms
    if(j==0) {int nHistos=num;}
    else {if(num != nHistos) {cerr << "ERROR: Input Root Files DO NOT have the same histograms!!" << endl;exit(1);} }
  }

  /* it's possible that not all root files were created by running over the same number of events. Therefore, the code below takes
     care of this by renomarlizing everything - renormalized events = max_events_analyzed * cut_efficiency. To do this, the maximum
     number of events analyzed must be determined. */
  float MaxEventsAnalyzed = 0;
  for(int i=0; i<theEventsAnalyzed.size(); i++) {
    if(theEventsAnalyzed.at(i) > MaxEventsAnalyzed) {MaxEventsAnalyzed = theEventsAnalyzed.at(i);}
  }

  //  outFile << "Detailed Efficiencies " << "\n";
  outFile << "Table of Efficiencies " << "\n\n";
  outFile << "Efficiencies below are calculated with respect to the number of events after skimming" << "\n\n";
  outFile << "-------------------------------------------------------------------------------------------------------------------------------\n";
  outFile << "         Name                         Events               Relative (%)                     Cumulative (%)\n";
  outFile << "-------------------------------------------------------------------------------------------------------------------------------\n";

  // calculate relative and cumulative cut efficiencies
  float effy = 1;
  float deffyOvereffy = 0;
  for (int i=0;i<theCumulativeEfficiency.size(); i++) {
    float numerator = (int)((MaxEventsAnalyzed * theCumulativeEfficiency.at(i))+0.5);
    if(i==0) {float denominator = (int)(MaxEventsAnalyzed + 0.5);}
    else {float denominator = (int)((MaxEventsAnalyzed * theCumulativeEfficiency.at(i-1))+0.5);}
//    float theRelativeEfficiency = numerator / denominator;
//    float efferror = sqrt(theRelativeEfficiency*(1.-theRelativeEfficiency)/denominator); // binomial uncertainty
    float theRelativeEfficiency;
    float efferror;
    if((numerator > 0) && (denominator > 0)) {
      theRelativeEfficiency = numerator / denominator;
      efferror = sqrt(theRelativeEfficiency*(1.-theRelativeEfficiency)/denominator);
    } else if((numerator == 0) && (denominator > 0)) {
      theRelativeEfficiency = 0.0;
      efferror = 1.0 / denominator;
    } else if((numerator == 0) && (denominator == 0)) {
      theRelativeEfficiency = 0.0;
      efferror = 0.0;
    } else {
      theRelativeEfficiency = numerator / denominator;
      efferror = sqrt(theRelativeEfficiency*(1.-theRelativeEfficiency)/denominator);
    }
    /* binomial uncertainties do not work when the efficiency gets close to 1 or 0 (e.g. efficiency cannot 
       be 99 +- 2 because the efficiency cannot be e.g. 101) ... in these cases, use bayesian */
    if(((theRelativeEfficiency + efferror) > 1.) || ((theRelativeEfficiency - efferror) < 0.)){
      TH1F* theNumHisto = new TH1F("theNumHisto","theNumHisto",1,0,1);
      theNumHisto->SetBinContent(1,numerator);
      theNumHisto->Sumw2();
      TH1F* theDenHisto = new TH1F("theDenHisto","",1,0,1);
      theDenHisto->SetBinContent(1,denominator);
      theDenHisto->Sumw2();
      TGraphAsymmErrors* bayesEff = new TGraphAsymmErrors();
      bayesEff->BayesDivide(theNumHisto,theDenHisto,"b");
      if(bayesEff->GetErrorYhigh(0) > bayesEff->GetErrorYlow(0)) {efferror = bayesEff->GetErrorYhigh(0);}
      else {efferror = bayesEff->GetErrorYlow(0);}
      delete theNumHisto;
      delete theDenHisto;
      delete bayesEff;
    }
    if(theRelativeEfficiency == 1.0) {efferror = 0;}
    // recalculate efficiencies and errors incorporating scale factors
    if((numerator > 0) && (denominator > 0)) {
      efferror = sqrt(pow(efferror/theRelativeEfficiency,2.0) + pow((inScaleFactorError.at(i)/inScaleFactor.at(i)),2.0));
      theRelativeEfficiency = theRelativeEfficiency * inScaleFactor.at(i);
      efferror = efferror * theRelativeEfficiency;
    } else if((numerator == 0) && (denominator > 0)) {
      efferror = inScaleFactor.at(i) / denominator;
      theRelativeEfficiency = theRelativeEfficiency * inScaleFactor.at(i);
    } else if((numerator == 0) && (denominator == 0)) {
      efferror = 0.0;
      theRelativeEfficiency = theRelativeEfficiency * inScaleFactor.at(i);
    } else {
      efferror = sqrt(pow(efferror/theRelativeEfficiency,2.0) + pow((inScaleFactorError.at(i)/inScaleFactor.at(i)),2.0));
      theRelativeEfficiency = theRelativeEfficiency * inScaleFactor.at(i);
      efferror = efferror * theRelativeEfficiency;
    }
//    efferror = sqrt(pow(efferror/theRelativeEfficiency,2.0) + pow((inScaleFactorError.at(i)/inScaleFactor.at(i)),2.0));
//    theRelativeEfficiency = theRelativeEfficiency * inScaleFactor.at(i);
//    efferror = efferror * theRelativeEfficiency;
    numerator = (int)((MaxEventsAnalyzed * theCumulativeEfficiency.at(i))+0.5);
    denominator = (int)(MaxEventsAnalyzed + 0.5);
//    float cumulativeEfficiency = numerator / denominator;
//    float efferror2 = sqrt(cumulativeEfficiency*(1.-cumulativeEfficiency)/denominator);
    float cumulativeEfficiency;
    float efferror2;
    if((numerator > 0) && (denominator > 0)) {
      cumulativeEfficiency = numerator / denominator;
      efferror2 = sqrt(cumulativeEfficiency*(1.-cumulativeEfficiency)/denominator);
    } else if((numerator == 0) && (denominator > 0)) {
      cumulativeEfficiency = 0.0;
      efferror2 = 1.0 / denominator;
    } else if((numerator == 0) && (denominator == 0)) {
      cumulativeEfficiency = 0.0;
      efferror2 = 0.0;
    } else {
      cumulativeEfficiency = numerator / denominator;
      efferror2 = sqrt(cumulativeEfficiency*(1.-cumulativeEfficiency)/denominator);
    }
    /* binomial uncertainties do not work when the efficiency gets close to 1 or 0 (e.g. efficiency cannot 
       be 99 +- 2 because the efficiency cannot be e.g. 101) ... in these cases, use bayesian */
    if(((cumulativeEfficiency + efferror2) > 1.) || ((cumulativeEfficiency - efferror2) < 0.)){
      TH1F* theNumHisto = new TH1F("theNumHisto","theNumHisto",1,0,1);
      theNumHisto->SetBinContent(1,numerator);
      theNumHisto->Sumw2();
      TH1F* theDenHisto = new TH1F("theDenHisto","",1,0,1);
      theDenHisto->SetBinContent(1,denominator);
      theDenHisto->Sumw2();
      TGraphAsymmErrors* bayesEff = new TGraphAsymmErrors();
      bayesEff->BayesDivide(theNumHisto,theDenHisto,"b");
      if(bayesEff->GetErrorYhigh(0) > bayesEff->GetErrorYlow(0)) {efferror2 = bayesEff->GetErrorYhigh(0);}
      else {efferror2 = bayesEff->GetErrorYlow(0);}
      delete theNumHisto;
      delete theDenHisto;
      delete bayesEff;
    }
    if(cumulativeEfficiency == 1.0) {efferror2 = 0;}
    // recalculate efficiencies and errors incorporating scale factors
    if((numerator > 0) && (denominator > 0)) {
      for(int numberOfSFs = 0; numberOfSFs <= i; numberOfSFs++) {
        efferror2 = sqrt(pow(efferror2/cumulativeEfficiency,2.0) + pow((inScaleFactorError.at(numberOfSFs)/inScaleFactor.at(numberOfSFs)),2.0));
        cumulativeEfficiency = cumulativeEfficiency * inScaleFactor.at(numberOfSFs);
        efferror2 = efferror2 * cumulativeEfficiency;
      }
    } else if((numerator == 0) && (denominator > 0)) {
      efferror2 = 1.0;
      for(int numberOfSFs = 0; numberOfSFs <= i; numberOfSFs++) {
        efferror2 = efferror2 * inScaleFactor.at(numberOfSFs);
        cumulativeEfficiency = cumulativeEfficiency * inScaleFactor.at(numberOfSFs);
      }
      efferror2 = efferror2 / denominator;
    } else if((numerator == 0) && (denominator == 0)) {
      efferror2 = 0.0;
      cumulativeEfficiency = 0.0;
    } else {
      for(int numberOfSFs = 0; numberOfSFs <= i; numberOfSFs++) {
        efferror2 = sqrt(pow(efferror2/cumulativeEfficiency,2.0) + pow((inScaleFactorError.at(numberOfSFs)/inScaleFactor.at(numberOfSFs)),2.0));
        cumulativeEfficiency = cumulativeEfficiency * inScaleFactor.at(numberOfSFs);
        efferror2 = efferror2 * cumulativeEfficiency;
      }
    }
//    for(int numberOfSFs = 0; numberOfSFs <= i; numberOfSFs++) {
//      efferror2 = sqrt(pow(efferror2/cumulativeEfficiency,2.0) + pow((inScaleFactorError.at(numberOfSFs)/inScaleFactor.at(numberOfSFs)),2.0));
//      cumulativeEfficiency = cumulativeEfficiency * inScaleFactor.at(numberOfSFs);
//      efferror2 = efferror2 * cumulativeEfficiency;
//    }
//    efferror2 = sqrt(pow(efferror2/cumulativeEfficiency,2.0) + pow((inScaleFactorError.at(i)/inScaleFactor.at(i)),2.0));
//    cumulativeEfficiency = cumulativeEfficiency * inScaleFactor.at(i);
//    efferror2 = efferror2 * cumulativeEfficiency;
    // output for cut flow efficiency table - efficiencies and errors
    outFile	<<setw(24)<<inProcess.at(i)
        	<<setw(20)<<(int)((MaxEventsAnalyzed * theCumulativeEfficiency.at(i))+0.5)
        	<<setw(15)<<setprecision(4)<<theRelativeEfficiency*100.0<<setw(4)<<" +- "<<setw(10)<<setprecision(4)<<(efferror*100.0)
        	<<setw(20)<<setprecision(4)<<cumulativeEfficiency*100.0<<setw(4)<<" +- "<<setw(10)<<setprecision(4)<<(efferror2 * 100.0)
        	<<endl;
    RelativeEffVector.push_back(theRelativeEfficiency);
    RelativeEffErrorVector.push_back(efferror);
    effy = effy * theRelativeEfficiency;
    deffyOvereffy = sqrt(pow(efferror/theRelativeEfficiency,2.0) + pow(deffyOvereffy,2.0));
    // output for cut flow efficiency table - calculation of the expected # of events at specific luminosity
    if(i==(theCumulativeEfficiency.size() - 1)) {
      outFile << "-------------------------------------------------------------------------------------------------------------------------------\n";
      outFile	<<setw(24)<<"Expected # of Events"
        	<<setw(20)<<"@ L ="
                <<setw(21)<<lumi<<" ipb "
        	<<setw(15)<<setprecision(4)<<"---> "<<(x_section.at(i) * lumi * effSkim * effy)
//                <<setw(4)<<" +- "<<setw(6)<<setprecision(4)<<(x_section.at(i) * lumi * effSkim * effy * deffyOvereffy)
                <<setw(4)<<" +- "<<setw(6)<<setprecision(4)<<(x_section.at(i) * lumi * effSkim * efferror2)
        	<<endl;
    }
  }
  outFile << "-------------------------------------------------------------------------------------------------------------------------------\n";

  if(( nHistList % inRootFiles.size() == 0 )) {
    for(int i=0; i<nHistos; i++) {
      string histname = HistList->At(i)->GetName();
      string histoEffyName = "hEffy_" + inProcess.at(0);
      TH1F *h = (TH1F*)HistList->At(i);
      h->Sumw2();
      if(IsData == "0") {
        if(h->Integral(0,(h->GetXaxis()->GetNbins()+1)) > 0) {h->Scale(1.0 / theEventsPassing.at(0) * x_section.at(0) * lumi * effSkim);}
        hEffy = new TH1F(histoEffyName.c_str(),histoEffyName.c_str(),
                         h->GetXaxis()->GetNbins(),h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
        for(int ibin=0; ibin<=(hEffy->GetXaxis()->GetNbins() + 1); ibin++) {
          hEffy->SetBinContent(ibin,RelativeEffVector.at(0));
          hEffy->SetBinError(ibin,0);
        }
        h->Multiply(h,hEffy,1,1,"");
        for(int ibin=0; ibin<=(h->GetXaxis()->GetNbins() + 1); ibin++) {
          h->SetBinError(ibin, sqrt(pow(h->GetBinError(ibin),2.0) + h->GetBinContent(ibin)) ); // propagating MC statistical uncertainty & poisson uncertainty on the rate
        }
      }
      h->SetName(histname.c_str());
      string YtitleString = "N / " + lumi_string + " pb^{-1}";
      h->GetYaxis()->SetTitle(YtitleString.c_str());
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetTitleFont(62);
      h->GetYaxis()->CenterTitle();
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelFont(62);
      h->GetXaxis()->SetTitle(theHistNameStrings.at(i).c_str());
      h->GetXaxis()->SetTitleSize(0.06);
      h->GetXaxis()->SetTitleFont(62);
      h->GetXaxis()->CenterTitle();
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetLabelFont(62);
      TFile *hfile = (TFile*)gROOT->FindObject(outputRootFileForNormalizedHistos.c_str());
      if (hfile) {hfile->Close();}
      hfile = new TFile(outputRootFileForNormalizedHistos.c_str(),"UPDATE");
      h->Write();
      hfile->Close();
      if(IsData == "0") {delete hEffy;}
      delete h;
      if(inRootFiles.size() > 1) {
        for(int j=1;j<inRootFiles.size();j++) {
          string hist2name = HistList->At(i+(j*nHistos))->GetName();
          if(theHistNameStrings.at(i) == theHistNameStrings.at(i+(j*nHistos))) {
            string histoEffyName2 = "hhEffy_" + inProcess.at(j);
            TH1F *hh = (TH1F*)HistList->At(i+(j*nHistos));
            hh->Sumw2();
            if(IsData == "0") {
              if(hh->Integral(0,(hh->GetXaxis()->GetNbins()+1)) > 0) {hh->Scale(1.0 / theEventsPassing.at(j) * x_section.at(j) * lumi * effSkim);}
              for(int jfile=0;jfile<=j;jfile++) {
                hhEffy = new TH1F(histoEffyName2.c_str(),histoEffyName2.c_str(),
                                  hh->GetXaxis()->GetNbins(),hh->GetXaxis()->GetXmin(),hh->GetXaxis()->GetXmax());
                for(int ibin=0; ibin<=(hhEffy->GetXaxis()->GetNbins() + 1); ibin++) {
                  hhEffy->SetBinContent(ibin,RelativeEffVector.at(jfile));
                  hhEffy->SetBinError(ibin,0);
                }
                hh->Multiply(hh,hhEffy,1,1,"");
                delete hhEffy;
              }
              for(int ibin=0; ibin<=(hh->GetXaxis()->GetNbins() + 1); ibin++) {
                hh->SetBinError(ibin, sqrt(pow(hh->GetBinError(ibin),2.0) + hh->GetBinContent(ibin)) );
              }
            }
            hh->SetName(hist2name.c_str());
            string YtitleString = "N / " + lumi_string + " pb^{-1}";
            hh->GetYaxis()->SetTitle(YtitleString.c_str());
            hh->GetYaxis()->SetTitleSize(0.06);
            hh->GetYaxis()->SetTitleFont(62);
            hh->GetYaxis()->CenterTitle();
            hh->GetYaxis()->SetLabelSize(0.05);
            hh->GetYaxis()->SetLabelFont(62);
            hh->GetXaxis()->SetTitle(theHistNameStrings.at(i+(j*nHistos)).c_str());
            hh->GetXaxis()->SetTitleSize(0.06);
            hh->GetXaxis()->SetTitleFont(62);
            hh->GetXaxis()->CenterTitle();
            hh->GetXaxis()->SetLabelSize(0.05);
            hh->GetXaxis()->SetLabelFont(62);
            TFile *hfile = (TFile*)gROOT->FindObject(outputRootFileForNormalizedHistos.c_str());
            if (hfile) {hfile->Close();}
            hfile = new TFile(outputRootFileForNormalizedHistos.c_str(),"UPDATE");
            hh->Write();
            hfile->Close();
	    delete hh;
          }
        }
      }
    }
  }

  std::cout << "I'm here ..." << std::endl;

  if(( nHistList2 % inRootFiles.size() == 0 )) {
    for(int i=0; i<nHistos; i++) {
      string histname = HistList2->At(i)->GetName();
      TH1F *h = (TH1F*)HistList2->At(i);
      h->Sumw2();
      if(h->Integral(0,(h->GetXaxis()->GetNbins()+1)) > 0) {h->Scale(1.0 / h->Integral(0,(h->GetXaxis()->GetNbins()+1)));}
      h->SetName(histname.c_str());
      h->GetYaxis()->SetTitle("a.u.");
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetTitleFont(62);
      h->GetYaxis()->CenterTitle();
      h->GetYaxis()->SetLabelSize(0.05);
      h->GetYaxis()->SetLabelFont(62);
      h->GetXaxis()->SetTitle(theHistNameStrings.at(i).c_str());
      h->GetXaxis()->SetTitleSize(0.06);
      h->GetXaxis()->SetTitleFont(62);
      h->GetXaxis()->CenterTitle();
      h->GetXaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetLabelFont(62);
      TFile *hfile = (TFile*)gROOT->FindObject(outputRootFileForProbabilityHistos.c_str());
      if (hfile) {hfile->Close();}
      hfile = new TFile(outputRootFileForProbabilityHistos.c_str(),"UPDATE");
      h->Write();
      hfile->Close();
      delete h;
      if(inRootFiles.size() > 1) {
        for(int j=1;j<inRootFiles.size();j++) {
          string hist2name = HistList2->At(i+(j*nHistos))->GetName();
          if(theHistNameStrings.at(i) == theHistNameStrings.at(i+(j*nHistos))) {
            TH1F *hh = (TH1F*)HistList2->At(i+(j*nHistos));
            hh->Sumw2();
            if(hh->Integral(0,(hh->GetXaxis()->GetNbins()+1)) > 0) {hh->Scale(1.0 / hh->Integral(0,(hh->GetXaxis()->GetNbins()+1)));}
            hh->SetName(hist2name.c_str());
            hh->GetYaxis()->SetTitle("a.u.");
            hh->GetYaxis()->SetTitleSize(0.06);
            hh->GetYaxis()->SetTitleFont(62);
            hh->GetYaxis()->CenterTitle();
            hh->GetYaxis()->SetLabelSize(0.05);
            hh->GetYaxis()->SetLabelFont(62);
            hh->GetXaxis()->SetTitle(theHistNameStrings.at(i+(j*nHistos)).c_str());
            hh->GetXaxis()->SetTitleSize(0.06);
            hh->GetXaxis()->SetTitleFont(62);
            hh->GetXaxis()->CenterTitle();
            hh->GetXaxis()->SetLabelSize(0.05);
            hh->GetXaxis()->SetLabelFont(62);
            TFile *hfile = (TFile*)gROOT->FindObject(outputRootFileForProbabilityHistos.c_str());
            if (hfile) {hfile->Close();}
            hfile = new TFile(outputRootFileForProbabilityHistos.c_str(),"UPDATE");
            hh->Write();
            hfile->Close();
	    delete hh;
          }
        }
      }
    }
  }

}





