{

  #include <iostream>
  #include <iomanip>
  #include <vector>
  gROOT->Reset();
  gStyle->SetOptTitle(0);

  // determine the current directory path
  string currentDirectory = gSystem->pwd();

  // open input file containing all information: folder names
  ifstream inFileFile;
  char inputFileFilename[] = "main.in";
  inFileFile.open(inputFileFilename, ios::in);
  // if can't open input file, exit the code
  if (!inFileFile) {
    cerr << "Can't open input file " << inputFileFilename << endl;
    exit(1);
  }

  vector<string> inputFolders;
  inputFolders.clear();
  string inputStringString;
  string inputTypeType;

  // grab all relevant information from the input file
  while (inFileFile >> inputTypeType >> inputStringString) {
    if(inputTypeType == "directory") {
      inputFolders.push_back(inputStringString); // directory path for each BG
    } else {
      cerr << "Incorrect input type " << inputTypeType << endl; // exit code if unwanted input is specified
      exit(1);
    }
  }

  std::cout << inputFolders.size() << std::endl;

  // loop over directories
  for(int jfolder=0;jfolder<inputFolders.size();jfolder++) {
    std::cout << jfolder << std::endl;
    gSystem->cd(inputFolders.at(jfolder).c_str());
    std::cout << gSystem->pwd() << std::endl;
    gROOT->ProcessLine(".x Normalizer.C");
    gSystem->cd(currentDirectory.c_str());
    std::cout << gSystem->pwd() << std::endl;
    std::cout << jfolder << std::endl;
  }

  gROOT->ProcessLine(".x Plotter.C");

}
