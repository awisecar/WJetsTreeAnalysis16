{
    string srcdir = "SourceFiles/";
    vector<string> sources;
    sources.push_back("getFilesAndHistograms");
    sources.push_back("functions");
    sources.push_back("Plotter");
    // sources.push_back("Plotter_ONLY_DATA");
    //--- Load shared libraries ---
    unsigned int nSources = sources.size();
    for (unsigned int i(0); i < nSources; i++){
        cout <<"Compiling " << srcdir + sources[i] << ".cc" << endl;
        gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc+").c_str());
    }

    // dit bonjour
    welcomeMessage();

    // plot only data -----
    // Plotter_ONLY_DATA("SMu", 2017, 30, 0, 0, 0, 0, -1, 0, -999999, 999999, 0, 0, 1); //no MET cut

    // plot data and MC -----
    // Plotter("SMu", 2017, 30, 0, 0, 0, 0, -1, 0, -999999, 999999, 0, 0, 1); //no MET cut
    // //Plotter("SMu", 2017, 30, 0, 0, 0, 30, -1, 0, -999999, 999999, 0, 0, 1); // METcut of 30 GeV

    // Plotter("SMu", 2017, 30, 0, 0, 0, 0, -1, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = -1
    // Plotter("SMu", 2017, 30, 0, 0, 0, 0,  2, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = 2
    Plotter("SMu", 2017, 30, 0, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements
    
    // grab statistics for the plotted data -----
    // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, -1, false); //no MET cut, doBJets = -1, doTTScale = False
    getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 0, false); //no MET cut, doBJets=0, doTTScale=False
    // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 2, false); //no MET cut, doBJets=2, doTTScale=False
    // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 30, -1, false); //30 GeV MET cut, doBJets = -1, doTTScale = False
}
