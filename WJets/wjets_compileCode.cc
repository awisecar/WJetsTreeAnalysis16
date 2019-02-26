{
    string srcdir = "Sources/";

//    //--- clean the old *_cc.d and *_cc.so files to allow for fresh compilation ---
//    //the "system" lines allow the script to run shell commands without using shell directly
//    string cmd = "if ls *_cc.d &> .ls_tmp.list; then rm *_cc.d; fi";
//    system(cmd.c_str());
//    cmd = "if ls *_cc.so &> .ls_tmp.list; then rm *_cc.so; fi";
//    system(cmd.c_str());
//    cmd = "if ls " + srcdir + "*_cc.d &> .ls_tmp.list; then rm " + srcdir + "*_cc.d; fi";
//    system(cmd.c_str());
//    cmd = "if ls " + srcdir + "*_cc.so &> .ls_tmp.list; then rm " + srcdir + "*_cc.so; fi";
//    system(cmd.c_str());
//    system("rm .ls_tmp.list");


    //--- compile the codes again
    vector<string> sources;
    sources.push_back("getFilesAndHistograms");
    sources.push_back("functions");
    sources.push_back("funcReweightResp");
    sources.push_back("HistoSet");
    sources.push_back("ZJetsAndDPS");
    unsigned int nSources = sources.size();
    gROOT->ProcessLine(".L /cvmfs/cms.cern.ch/slc5_amd64_gcc434/external/lhapdf/5.8.5/lib/libLHAPDF.so");
    for (unsigned int i(0); i < nSources; i++) {
        std::cout << "Compiling " << srcdir + sources[i] << ".cc" << std::endl;
        //need the plus at the end of ".cc" in this line since we're using CINT in this version of ROOT
	gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc+").c_str());
    }

}
