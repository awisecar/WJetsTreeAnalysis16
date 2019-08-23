void rootlogon(void) {
  cout << "\nExecuting rootlogon.C" << endl;
  string currentFile =  __FILE__;
  string currentWorkingDir = currentFile.substr(0, currentFile.find("./rootlogon.C"));
  cout << "Current working directory is:\n\t" << currentWorkingDir << endl;

  // gErrorIgnoreLevel = kError;
  
  // cout << "Adding " << incdir << " to includes directories..." << endl;
  // gSystem->AddIncludePath(string(" -I" + incdir + " ").c_str());

  string srcDir = "SourceFiles/";
  vector<string> sources;
  sources.push_back("getFilesAndHistograms");
  sources.push_back("functions");
  sources.push_back("HistoSet");
  sources.push_back("ZJetsAndDPS");
  unsigned int nSources = sources.size();
  for (unsigned int i(0); i < nSources; i++) {
    std::cout << " >>> Compiling/loading libraries for " << srcDir + sources[i] << ".cc" << std::endl;

    // single-plus at the end rebuilds the library only if the script 
    // or any of the files it includes are newer than the library
    // gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc+").c_str());
    gROOT->LoadMacro(string(srcDir + sources[i] + ".cc+").c_str());
  }
}
