{
    string srcdir = "SourceFiles/";
    vector<string> sources;
    sources.push_back("getFilesAndHistograms");
    sources.push_back("functions");
    sources.push_back("Plotter");
    //--- Load shared libraries ---
    unsigned int nSources = sources.size();
    for (unsigned int i(0); i < nSources; i++){
        cout <<"Compiling " << srcdir + sources[i] << ".cc" << endl;
        gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc+").c_str());
    }

    Plotter("SMu", 30, 0, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut
   
    //last option of function is doTTScale
    getStatistics("SMu", 30, 0, false, true, 0, false, false, 0, 0, false); //no MET cut, doBJets = 0, doTTScale = False

}
