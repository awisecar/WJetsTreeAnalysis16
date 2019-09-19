{
    string srcdir = "SourceFiles/";
    vector<string> sources;
    sources.push_back("getFilesAndHistograms");
    sources.push_back("functions");
    sources.push_back("DataDrivenQCD");
    //--- Load shared libraries ---
    unsigned int nSources = sources.size();
    for (unsigned int i(0); i < nSources; i++){
        cout <<"Compiling " << srcdir + sources[i] << ".cc" << endl;
        gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc+").c_str());
    }

    //arguments are leptonFlavor, METcut, doBJets
    // DataDrivenQCD("SMu", 0, -1); // bveto on 1 jet
    DataDrivenQCD("SMu", 0, 0); // no bveto

    //DataDrivenQCD("SMu", 30, -1); // METcut of 30 GeV

}
