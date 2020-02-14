void runDYJets(int doWhat = 0, int doQCD = 0, int doSysRunning = 0)
{
    string srcdir = "SourceFiles/";
    vector<string> sources;
    sources.push_back("getFilesAndHistograms");
    sources.push_back("functions");
    sources.push_back("HistoSet");
    sources.push_back("ZJetsAndDPS");
    ////--- Load shared libraries (the .so's) ---
    unsigned int nSources = sources.size();
    for (unsigned int i(0); i < nSources; i++) {
        std::cout << "Compiling " << srcdir + sources[i] << ".cc" << std::endl;
        gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc+").c_str());
    }
        
    //------------------
     // int doWhat       = 41;
                              // 100 - all ; 10, 11, ... - individual data samples, 1 - background , 2 - tau ?, 3 - DY, 
                              // 41 - W+jets inc. NLO-FxFx, 42 - W+jets inc. LO-MLM
                              // 5 - W+jets FxFx W pT-binned, 6 - W+jets FxFx jet-binned,
                              // 51 - MC gen, 90 - PDF Syst., 1001 - do pull DY samples
        
     // int doSysRunning = 0;
                             // 0 - no syst running, 100 - all systematic runnings,
                             // 1 - PU, 2 - JES, 3 - XSEC, 4 - JER, 5 - LepSF,
                             // 6 - BtagSF, 7 - MES, 8 - MER, 9 - WB, 10 - RESP
        
     // int doQCD        = 0;
                             // 0-3 : 4 combination between isolation/anti-isolation and MT cuts for QCD BG estimation
        
    int doBJets      = -1; //normal btag veto
    //int doBJets      = 2; //ttbar SFs
                            // 0 - no infor on B-jets will be used ;
                            // 1, 2 .. require at least 1, 2, .. ; use 2 for ttbar systmatics;
                            // -1, -2, .. veto the event if you have 1 or more, 2 or more .. b-jets ;
                            // 101 - require exactly 1 b-jet
    //------------------
    
    string lepSelection = "SMu";
    double muLumi(35916.625); // 80X 2016 data bonzai 23Sep2016ReReco golden json

    //switches
    bool doSSign  =  0;     // contribution of QCD to emu in TTbar
    bool doInvMassCut = 0;
    int doPUStudy = -10 ;  // default int the ZJets
    bool doFlat   = 0;
    bool doRoch   = 0;
    bool doVarWidth = 1;
    bool hasPartonInfo = 0;
    bool doDataEff(0);

    //phase space cuts
    int jetPtMin = 30;
    int jetPtMax = 0;      // 0 - for no jetPtMax cut
    int ZPtMin = 0;
    int ZEtaMin  = -999999;  // default value -999999  !!!!!!!  factor 100 to keep things integer .... eta 2.4  = eta Cut 240
    int ZEtaMax  = 999999;   // default value  999999
    int METcut = 0; 
    int jetEtaMin = -24;
    int jetEtaMax = 24;

    const int NSystData(3), NSystMC(12); // set NSystMC to 5 if you want to do only PU, XSEC
    const int NSystWJets(16), NSystDYJets(14);

    //systematic variation tags for data    
    short dataSyst[NSystData]  = {0, 2, 2};
    short dataDir[NSystData]   = {0,-1, 1};

    //systematic variation tags for background (generic)
    short bgSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short bgDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float bgScale[NSystMC]     = {1, 1, 1,   0.,   0.,  1, 1, 1, 1, 1, 1, 1};
 
    //systematic variations for background (ttbar, diboson, single top, DY+jets)
    short ttSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short ttDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float ttScale[NSystMC]     = {1, 1, 1,  0.06,  0.06,  1, 1, 1, 1, 1, 1, 1};
    
    short zzSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short zzDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float zzScale[NSystMC]     = {1, 1, 1, 0.07, 0.07,  1, 1, 1, 1, 1, 1, 1};
        
    short wzSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short wzDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float wzScale[NSystMC]     = {1, 1, 1, 0.07, 0.07,  1, 1, 1, 1, 1, 1, 1};
    
    short wwSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short wwDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float wwScale[NSystMC]     = {1, 1, 1, 0.06, 0.06,  1, 1, 1, 1, 1, 1, 1};
        
    short tcsSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short tcsDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float tcsScale[NSystMC]    = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1};
    
    short tctSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short tctDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float tctScale[NSystMC]    = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1};
    
    short tcwSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short tcwDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float tcwScale[NSystMC]    = {1, 1, 1, 0.06, 0.06,  1, 1, 1, 1, 1, 1, 1};
        
    short dySyst[NSystDYJets]  = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 4, 4, 8};
    short dyDir[NSystDYJets]   = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1,-1, 1, 1};
    float dyScale[NSystDYJets] = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1};

    //systematic variations for signal (W+jets)
    short wjSyst[NSystWJets]   = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 4, 4, 8, 9, 10};
    short wjDir[NSystWJets]    = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1,-1, 1, 1, 1,  1};
    float wjScale[NSystWJets]  = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1};
        
    //------------------

    // Data
    if (doWhat == 10 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

            ZJetsAndDPS DMudata1(lepSelection+"_13TeV_Data_dR_5311_List_1", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata1.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata2(lepSelection+"_13TeV_Data_dR_5311_List_2", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata2.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata3(lepSelection+"_13TeV_Data_dR_5311_List_3", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata3.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 11 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata4(lepSelection+"_13TeV_Data_dR_5311_List_4", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata4.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata5(lepSelection+"_13TeV_Data_dR_5311_List_5", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata5.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata6(lepSelection+"_13TeV_Data_dR_5311_List_6", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata6.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 12 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata7(lepSelection+"_13TeV_Data_dR_5311_List_7", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata7.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata8(lepSelection+"_13TeV_Data_dR_5311_List_8", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata8.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata9(lepSelection+"_13TeV_Data_dR_5311_List_9", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata9.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 13 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata10(lepSelection+"_13TeV_Data_dR_5311_List_10", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata10.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata11(lepSelection+"_13TeV_Data_dR_5311_List_11", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata11.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata12(lepSelection+"_13TeV_Data_dR_5311_List_12", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata12.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 14 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata13(lepSelection+"_13TeV_Data_dR_5311_List_13", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata13.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata14(lepSelection+"_13TeV_Data_dR_5311_List_14", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata14.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata15(lepSelection+"_13TeV_Data_dR_5311_List_15", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata15.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 15 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata16(lepSelection+"_13TeV_Data_dR_5311_List_16", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata16.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata17(lepSelection+"_13TeV_Data_dR_5311_List_17", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata17.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata18(lepSelection+"_13TeV_Data_dR_5311_List_18", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata18.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 16 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata19(lepSelection+"_13TeV_Data_dR_5311_List_19", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata19.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata20(lepSelection+"_13TeV_Data_dR_5311_List_20", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata20.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata21(lepSelection+"_13TeV_Data_dR_5311_List_21", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata21.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 17 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata22(lepSelection+"_13TeV_Data_dR_5311_List_22", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata22.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata23(lepSelection+"_13TeV_Data_dR_5311_List_23", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata23.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata24(lepSelection+"_13TeV_Data_dR_5311_List_24", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata24.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 18 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata25(lepSelection+"_13TeV_Data_dR_5311_List_25", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata25.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata26(lepSelection+"_13TeV_Data_dR_5311_List_26", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata26.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata27(lepSelection+"_13TeV_Data_dR_5311_List_27", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata27.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    if (doWhat == 19 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata28(lepSelection+"_13TeV_Data_dR_5311_List_28", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata28.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata29(lepSelection+"_13TeV_Data_dR_5311_List_29", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata29.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMudata30(lepSelection+"_13TeV_Data_dR_5311_List_30", 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMudata30.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }
        
    // Background -- TTbar
    if (doWhat == 21 || doWhat == 100 ){
        for (unsigned int i(0); i < NSystMC; i++){
            if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuTT(lepSelection+"_13TeV_TTJets_dR_5311_List", muLumi * 831.7  , 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuTT.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    // Background -- Diboson
    if (doWhat == 22 || doWhat == 100 ){
        for (unsigned int i(0); i < NSystMC; i++){
            if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuZZInc(lepSelection+"_13TeV_ZZ_dR_5311_List", muLumi * 15.4 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuZZInc.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMuWZInc(lepSelection+"_13TeV_WZ_dR_5311_List", muLumi * 23.5 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWZInc.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMuWWInc(lepSelection+"_13TeV_WW_dR_5311_List", muLumi * 12.21  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWWInc.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    // Background -- ST_s,ST_t
    if (doWhat == 23 || doWhat == 100 ){
        for (unsigned int i(0); i < NSystMC; i++){
            if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuT1(lepSelection+"_13TeV_ST_s_channel_dR_5311_List", muLumi * 3.35  ,  1., 1, !doDataEff, tcsSyst[i], tcsDir[i], tcsScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuT1.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMuT2(lepSelection+"_13TeV_ST_t_top_channel_dR_5311_List", muLumi * 136.02  ,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuT2.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMuT3(lepSelection+"_13TeV_ST_t_antitop_channel_dR_5311_List", muLumi * 80.95  ,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuT3.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    // Background -- ST_tW
    if (doWhat == 24 || doWhat == 100 ){
        for (unsigned int i(0); i < NSystMC; i++){
            if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuT4(lepSelection+"_13TeV_ST_tW_top_channel_dR_5311_List", muLumi * 35.85  , 1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuT4.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);

            ZJetsAndDPS DMuT5(lepSelection+"_13TeV_ST_tW_antitop_channel_dR_5311_List", muLumi * 35.85 ,  1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuT5.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    
    // DYJets_MIX_UNFOLDING_dR_5311_Inf3
    if (doWhat == 30 || doWhat == 100 ){
        int doGen = 0 ;
        if ( lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0 ) doGen = 1;
        
        for (unsigned int i(0); i < NSystDYJets; i++){
            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && dySyst[i] == 4) continue; // jet smearing part -- not done for SMu ---
            if ((lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0) && dySyst[i] == 3) continue; // xsec -- not done for Z+jets ---
            if (dySyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuDYMix(lepSelection+"_13TeV_DYJets50toInf_dR_5311_List", muLumi * 5765.4 , 1., 1, !doDataEff, dySyst[i], dyDir[i], dyScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuDYMix.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
        
    }
 
    //W+jets inclusive NLO-FxFx sample   
    if (doWhat == 41 || doWhat == 100){
        int doGen = 0;
        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0) && lepSelection.find("SMuE") == -1)  doGen = 1 ;

        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuWJFxFx(lepSelection+"_13TeV_WJets_FxFx_dR_5311_List", muLumi* 60290.0 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
        }
    }

    //W+jets inclusive LO-MLM sample
    if (doWhat == 42 || doWhat == 100 ){
       int doGen = 0;
       if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;

       for (unsigned int i(0); i < NSystWJets; i++){
           if (!doGen ) continue;
           if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
           if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

           ZJetsAndDPS DMuWJMLM(lepSelection+"_13TeV_WJets_MLM_dR_5311_List", muLumi* 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
           DMuWJMLM.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo);
       }
    }    
    
    // W+jets FxFx W pT-binned signal sample - 0 to 50 pT
    if (doWhat == 51 || doWhat == 100 ){
        int doGen = 0 ;
        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuWJFxFx_Wpt1(lepSelection+"_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List", muLumi* 56306.4 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_Wpt1.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization
        }
    }

    // W+jets FxFx W pT-binned signal sample - 50 to 100 pT
    if (doWhat == 52 || doWhat == 100 ){
        int doGen = 0 ;
        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuWJFxFx_Wpt2(lepSelection+"_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List", muLumi* 3241.33 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_Wpt2.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization
        }
    }

    // W+jets FxFx W pT-binned signal sample - 100 to 250 pT
    if (doWhat == 53 || doWhat == 100 ){
        int doGen = 0 ;
        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuWJFxFx_Wpt3(lepSelection+"_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List", muLumi* 677.82 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_Wpt3.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization
        }
    }

    // W+jets FxFx W pT-binned signal sample - 250 to Inf pT
    if (doWhat == 54 || doWhat == 100 ){
        int doGen = 0 ;
        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuWJFxFx_Wpt4(lepSelection+"_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List", muLumi* 24.083 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_Wpt4.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization

            ZJetsAndDPS DMuWJFxFx_Wpt5(lepSelection+"_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List", muLumi* 3.0563 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_Wpt5.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization
            
            ZJetsAndDPS DMuWJFxFx_Wpt6(lepSelection+"_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List", muLumi* 0.4602 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_Wpt6.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization
        }
    }

    
    // W+jets FxFx jet-binned signal sample - W+0J
    if ( doWhat == 61 || doWhat == 100 ){
        int doGen = 0 ;
        if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", muLumi* 49264.92 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_jet0.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization
        }
    }

    // W+jets FxFx jet-binned signal sample - W+1J
    if ( doWhat == 62 || doWhat == 100 ){
        int doGen = 0 ;
        if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", muLumi* 8280.36 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_jet1.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization

        }
    }

    // W+jets FxFx jet-binned signal sample - W+2J
    if ( doWhat == 63 || doWhat == 100 ){
        int doGen = 0 ;
        if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", muLumi* 3118.08 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, 0, jetEtaMin, jetEtaMax);
            DMuWJFxFx_jet2.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth, hasPartonInfo); //FxFx with NLO normalization
        }
    }

}
