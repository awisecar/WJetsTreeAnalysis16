{
    string srcdir = "Sources/";

    vector<string> sources;
    sources.push_back("getFilesAndHistograms");
    sources.push_back("functions");
    sources.push_back("Plotter");

    //--- Load shared libraries ---
    unsigned int nSources = sources.size();
    for (unsigned int i(0); i < nSources; i++){
        cout <<"Compiling " << srcdir + sources[i] << ".cc" << endl;
        gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc++").c_str());
    }

    //Plotter("SMu",30);
    // Plotter("SMu", 30, 0, 0, 0, 0, -1, 0, -999999, 999999, 0, 0, 1);
    Plotter("SMu", 30, 0, 0, 0, 30, -1, 0, -999999, 999999, 0, 0, 1); // METcut of 30 GeV
   
    //last option of function is doTTScale
    //getStatistics("SMu", 30, 0, false, true, 0, false, false, 0, -1, true);
    // getStatistics("SMu", 30, 0, false, true, 0, false, false, 0, -1, false);

    //getStatistics("SMu", 30, 0, false, true, 0, false, false, 0, 0);

    //--- clean the *_cc.d and *_cc.so files ---
    string cmd = "if ls *_cc.d &> .ls_tmp.list; then rm *_cc.d; fi";
    system(cmd.c_str());
    cmd = "if ls *_cc.so &> .ls_tmp.list; then rm *_cc.so; fi";
    system(cmd.c_str());
    cmd = "if ls " + srcdir + "*_cc.d &> .ls_tmp.list; then rm " + srcdir + "*_cc.d; fi";
    system(cmd.c_str());
    cmd = "if ls " + srcdir + "*_cc.so &> .ls_tmp.list; then rm " + srcdir + "*_cc.so; fi";
    system(cmd.c_str());
    system("rm .ls_tmp.list");

}