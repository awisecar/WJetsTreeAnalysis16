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

    //////////////////////////////////
    //         2016 data/MC         //
    //////////////////////////////////

    // Plotter("SMu", 2016, 30, 0, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=0
    // //Plotter("SMu", 2016, 30, 1, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=1
    // //Plotter("SMu", 2016, 30, 2, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=2
    // //Plotter("SMu", 2016, 30, 3, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=3

    // getStatistics("SMu", 2016, 30, 0, false, true, 0, false, false, 0, 0, false, true); //no MET cut, no b-tag requirements, doQCD=0, doTTScale=False, incl QCD BG
    // // getStatistics("SMu", 2016, 30, 0, false, true, 0, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=0, doTTScale=False, no QCD BG
    // // getStatistics("SMu", 2016, 30, 0, false, true, 1, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=1, doTTScale=False, no QCD BG
    // // getStatistics("SMu", 2016, 30, 0, false, true, 2, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=2, doTTScale=False, no QCD BG
    // // getStatistics("SMu", 2016, 30, 0, false, true, 3, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=3, doTTScale=False, no QCD BG

    // // for ttbar studies
    // before rescaling...
    // Plotter("SMu", 2016, 30, 0, 0, 0, 0,  2, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = 2
    // getStatistics("SMu", 2016, 30, 0, false, true, 0, false, false, 0, 2, false, false); // no MET cut, doBJets=2, doQCD=0, doTTScale=False, no QCD BG
    // after rescaling...
    // Plotter("SMu", 2016, 30, 0, 0, 0, 0,  2, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = 2
    // getStatistics("SMu", 2016, 30, 0, false, true, 0, false, false, 0, 2, true, false); // no MET cut, doBJets=2, doQCD=0, doTTScale=True, no QCD BG

    // running signal region with no btag veto after ttBar SFs applied
 //   Plotter("SMu", 2016, 30, 0, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=0
 //   getStatistics("SMu", 2016, 30, 0, false, true, 0, false, false, 0, 0, true, true); //no MET cut, no b-tag requirements, doQCD=0, doTTScale=True, incl QCD BG


    //////////////////////////////////
    //         2017 data/MC         //
    //////////////////////////////////

    // // plot only data -----
    // // Plotter_ONLY_DATA("SMu", 2017, 30, 0, 0, 0, 0, -1, 0, -999999, 999999, 0, 0, 1); //no MET cut

    // // plot data and MC -----
    // // //Plotter("SMu", 2017, 30, 0, 0, 0, 30, -1, 0, -999999, 999999, 0, 0, 1); // METcut of 30 GeV

    // // Plotter("SMu", 2017, 30, 0, 0, 0, 0, -1, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = -1
    // // Plotter("SMu", 2017, 30, 0, 0, 0, 0,  2, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = 2

    // Plotter("SMu", 2017, 30, 0, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=0
    // // Plotter("SMu", 2017, 30, 1, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=1
    // // Plotter("SMu", 2017, 30, 2, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=2
    // // Plotter("SMu", 2017, 30, 3, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=3
    
    // // grab statistics for the plotted data -----
    // // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, -1, false); //no MET cut, doBJets = -1, doTTScale = False
    // // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 0, false); //no MET cut, doBJets=0, doTTScale=False
    // // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 2, false); //no MET cut, doBJets=2, doTTScale=False
    // // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 30, -1, false); //30 GeV MET cut, doBJets = -1, doTTScale = False

    // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 0, false, true); //no MET cut, no b-tag requirements, doQCD=0, doTTScale=False, incl QCD BG
    // // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=0, doTTScale=False, no QCD BG
    // // getStatistics("SMu", 2017, 30, 0, false, true, 1, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=1, doTTScale=False, no QCD BG
    // // getStatistics("SMu", 2017, 30, 0, false, true, 2, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=2, doTTScale=False, no QCD BG
    // // getStatistics("SMu", 2017, 30, 0, false, true, 3, false, false, 0, 0, false, false); // no MET cut, no b-tag requirements, doQCD=3, doTTScale=False, no QCD BG


    // for ttbar studies
    // before rescaling...
    // Plotter("SMu", 2017, 30, 0, 0, 0, 0,  2, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = 2
    // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 2, false, false); // no MET cut, doBJets=2, doQCD=0, doTTScale=False, no QCD BG
    // after rescaling...
    // Plotter("SMu", 2017, 30, 0, 0, 0, 0,  2, 0, -999999, 999999, 0, 0, 1); //no MET cut, doBJets = 2
    // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 2, true, false); // no MET cut, doBJets=2, doQCD=0, doTTScale=True, no QCD BG


    // running signal region with no btag veto after ttBar SFs applied
    // Plotter("SMu", 2017, 30, 0, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=0
    // getStatistics("SMu", 2017, 30, 0, false, true, 0, false, false, 0, 0, true, true); //no MET cut, no b-tag requirements, doQCD=0, doTTScale=True, incl QCD BG


    //////////////////////////////////
    //         2018 data/MC         //
    //////////////////////////////////

    // running signal region with no btag veto
    Plotter("SMu", 2018, 30, 0, 0, 0, 0, 0, 0, -999999, 999999, 0, 0, 1); //no MET cut, no b-tag requirements, doQCD=0
    getStatistics("SMu", 2018, 30, 0, false, true, 0, false, false, 0, 0, false, true); //no MET cut, no b-tag requirements, doQCD=0, doTTScale=False, incl QCD BG





}
