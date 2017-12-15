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

    //Note--> In Plotter, the function "getFile" is used to grab each of the necessary files (from directory FILESDIRECTORY, or "HistoFiles/"): the variable fileSelect is set equal to a vector of integers (e.g. FilesTTBarWJets) corresponding to the different entries of the structure ProcessInfo (in fileNames.h), which contains the names of the data and BG files we're using
    //If you want to change the files you're stacking the histograms of, do it in fileNames.h

    //Plotter("SMu",30);
    Plotter("SMu", 30, 0, 0, 0, 0, 2);

    //last option of function is doTTScale
    getStatistics("SMu", 30, 0, false, true, 0, false, false, 0, 2, true);
    //getStatistics("SMu", 30, 0, false, true, 0, false, false, 0, 2, false);

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
