{
  cout << "\nExecuting rootlogon.C" << endl;
  string currentFile =  __FILE__;
  string currentWorkingDir = currentFile.substr(0, currentFile.find("./rootlogon.C"));
  cout << "Current working directory is:\n\t" << currentWorkingDir << endl;
  gErrorIgnoreLevel = kError;
  string incdir = currentWorkingDir + "Includes/";
  cout << "Adding " << incdir << " to includes directories..." << endl;
  gSystem->AddIncludePath(string("-I" + incdir).c_str());
  //cout << "Include Path -D__USE_XOPEN2K8 to fix lxplus6 compatibility" << endl;
  //gSystem->AddIncludePath("-D__USE_XOPEN2K8");
  cout << "\n";
}
