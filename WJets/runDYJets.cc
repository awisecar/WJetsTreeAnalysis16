void runDYJets(int doWhat = 0, int doQCD = 0, int doSysRunning = 0)
{
    string srcdir = "Sources/";

    vector<string> sources;
    sources.push_back("getFilesAndHistograms");
    sources.push_back("functions");
    sources.push_back("funcReweightResp");
    sources.push_back("HistoSet");
    sources.push_back("ZJetsAndDPS");

    ////--- Load shared libraries (the .so's) ---
    unsigned int nSources = sources.size();
    gSystem->AddIncludePath("-D__USE_XOPEN2K8");
    //gROOT->ProcessLine(".L /usr/local/lib/libLHAPDF.dylib");
    gROOT->ProcessLine(".L /cvmfs/cms.cern.ch/slc5_amd64_gcc434/external/lhapdf/5.8.5/lib/libLHAPDF.so");
    for (unsigned int i(0); i < nSources; i++) {
        std::cout << "Compiling " << srcdir + sources[i] << ".cc" << std::endl;
        gROOT->ProcessLine(string(".L " + srcdir + sources[i] + ".cc+").c_str());
    }
        
    //------
    //int doWhat       = 100;
                              // 100 - all ; 10, 11, ... - individual data samples, 1 - background , 2 - tau ?, 3 - DY, 
                              // 41 - W+jets inc. NLO-FxFx, 42 - W+jets inc. LO-MLM
                              // 5 - W+jets FxFx W pT-binned, 6 - W+jets FxFx jet-binned,
                              // 51 - MC gen, 90 - PDF Syst., 1001 - do pull DY samples
        
    //int doSysRunning = 0;
                             // 0 - no syst running, 100 - all systematic runnings,
                             // 1 - PU, 2 - JES, 3 - XSEC, 4 - JER, 5 - LepSF,
                             // 6 - BtagSF, 7 - MES, 8 - MER, 9 - WB, 10 - RESP
        
    //int doQCD        = 0;
                             // 0-3 : 4 combination between isolation/anti-isolation and MT cuts for QCD BG estimation
        
    int doBJets      = -1;
                            // 0 - no infor on B-jets will be used ;
                            // 1, 2 .. require at least 1, 2, .. ; use 2 for ttbar systmatics;
                            // -1, -2, .. veto the event if you have 1 or more, 2 or more .. b-jets ;
                            // 101 - require exactly 1 b-jet
    
        
    //--- Internal configuration (no need to change this)---
        
    /*string lepSelection = "SMu"; // default lumi is set for double muon dataset
    double muLumi(19.549); // DoubleMu number with pixelCalc
    double eleLumi(19.602); // DoubleEle number with pixelCalc
    if      (lepSelection == "DE")   muLumi = 19.602;
    else if (lepSelection == "SMuE") muLumi = 19.673;
    else if (lepSelection == "SMu")  muLumi = 19.244;
    else if (lepSelection == "SE")   muLumi = 19.174;
    */
    
    string lepSelection = "SMu"; // default lumi is set for double muon dataset
    double muLumi(35916.637); // 80X 2016 data bonzai 23Sep2016ReReco golden json

    double w_sum_WJets(3.73654e+12);
    double w_sum_TTbar(1.15005e+07);
    double w_sum_DYJets(4.48707e+11);
    double w_sum_ST_s(3.3187e+06);
    //double w_sum_ST_t();
    double w_sum_ST_tW_top(995600);
    double w_sum_ST_tW_antitop(988500);
    double w_sum_ZZ(996944);
    double w_sum_WW(1.9652e+06);
    double w_sum_WZ(1.08822e+08);


    int doRoch   = 0;
    int doFlat   = 0;
    int doPUStudy = -10 ;  // default int the ZJets
    bool doSSign  =  0;     // contribution of QCD to emu in TTbar
    bool doInvMassCut = 0 ;
    bool doDataEff(0);
        
    int doMETcut = 0 ;     // if you want to apply MET cut
    int jetPtMin = 30;
    int jetPtMax = 0;      // 0 - for no jetPtMax cut
    int ZPtMin = 0 ;
    int ZEtaMin  = -999999;  // default value -999999  !!!!!!!  factor 100 to keep things integer .... eta 2.4  = eta Cut 240
    int ZEtaMax  = 999999;   // default value  999999

    const int NSystData(3), NSystMC(12); // set NSystMC to 5 if you want to do only PU, XSEC
    const int NSystWJets(16), NSystDYJets(14);
    
    short dataSyst[NSystData]  = {0, 2, 2};
    short dataDir[NSystData]   = {0,-1, 1};

    short ttSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short ttDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float ttScale[NSystMC]     = {1, 1, 1,  0.06,  0.06,  1, 1, 1, 1, 1, 1, 1};

    short tauSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short tauDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float tauScale[NSystMC]    = {1, 1, 1,   0.,   0.,  1, 1, 1, 1, 1, 1, 1};

    short bgSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short bgDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float bgScale[NSystMC]     = {1, 1, 1,   0.,   0.,  1, 1, 1, 1, 1, 1, 1};
    
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
    float tctScale[NSystMC]    = {1, 1, 1, 0.05, 0.05,  1, 1, 1, 1, 1, 1, 1};
    
    short tcwSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8};
    short tcwDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1};
    float tcwScale[NSystMC]    = {1, 1, 1, 0.06, 0.06,  1, 1, 1, 1, 1, 1, 1};
        
    short wjSyst[NSystWJets]   = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 4, 4, 8, 9, 10};
    short wjDir[NSystWJets]    = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1,-1, 1, 1, 1,  1};
    float wjScale[NSystWJets]  = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1};
        
    short dySyst[NSystDYJets]  = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 4, 4, 8};
    short dyDir[NSystDYJets]   = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1,-1, 1, 1};
    float dyScale[NSystDYJets] = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1};
        

    // Data
    if (doWhat == 10 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;
                
            ZJetsAndDPS DMudata1(lepSelection+"_13TeV_Data_01_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata1.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 11 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata2(lepSelection+"_13TeV_Data_02_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata2.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 12 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata3(lepSelection+"_13TeV_Data_03_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata3.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 13 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata4(lepSelection+"_13TeV_Data_04_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata4.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 14 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata5(lepSelection+"_13TeV_Data_05_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata5.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 15 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata6(lepSelection+"_13TeV_Data_06_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata6.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 16 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata7(lepSelection+"_13TeV_Data_07_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata7.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 17 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata8(lepSelection+"_13TeV_Data_08_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata8.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 18 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata9(lepSelection+"_13TeV_Data_09_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata9.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }

    if (doWhat == 19 || doWhat == 100) {
        for (unsigned int i(0); i < NSystData; i++) {
            if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
            ZJetsAndDPS DMudata10(lepSelection+"_13TeV_Data_10_dR_5311_List", 1., 1, 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMudata10.Loop(1, 0, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy);
        }
    }
        
    // Background
    if (doWhat == 1 || doWhat == 100 ){
        for (unsigned int i(0); i < NSystMC; i++){
            if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuTT(lepSelection+"_13TeV_TTJets_dR_5311_List", muLumi * 831.7  , 1, 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin, ZEtaMax);
            DMuTT.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuZZInc(lepSelection+"_13TeV_ZZ_dR_5311_List", muLumi * 15.4 ,  1, 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,  ZEtaMax);
            DMuZZInc.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuWZInc(lepSelection+"_13TeV_WZ_dR_5311_List", muLumi * 23.5 , 1, 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin, ZEtaMax);
            DMuWZInc.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuWWInc(lepSelection+"_13TeV_WW_dR_5311_List", muLumi * 12.21  , 1, 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin, ZEtaMax);
            DMuWWInc.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuT1(lepSelection+"_13TeV_ST_s_channel_dR_5311_List", muLumi * 3.35  ,  1, 1, !doDataEff, tcsSyst[i], tcsDir[i], tcsScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin, ZEtaMax);
            DMuT1.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuT2(lepSelection+"_13TeV_ST_t_top_channel_dR_5311_List", muLumi * 136.02  ,  1, 1, !doDataEff, tcsSyst[i], tcsDir[i], tcsScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMuT2.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuT3(lepSelection+"_13TeV_ST_t_antitop_channel_dR_5311_List", muLumi * 80.95  ,  1, 1, !doDataEff, tcsSyst[i], tcsDir[i], tcsScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMuT3.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuT4(lepSelection+"_13TeV_ST_tW_top_channel_dR_5311_List", muLumi * 35.85  , 1, 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMuT4.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            ZJetsAndDPS DMuT5(lepSelection+"_13TeV_ST_tW_antitop_channel_dR_5311_List", muLumi * 35.85 ,  1, 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax);
            DMuT5.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);

            
        //    if ( (lepSelection.find("SE") == -1 && lepSelection.find("SMu") == -1 )  ) {
        //        // for Z plus jets: WJets is background
        //        ZJetsAndDPS DMuWJ(lepSelection+"_13TeV_WJetsALL_UNFOLDING_dR_5311",           muLumi*36703.   *1000/76102995.,1, 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuWJ.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //        
        //        ZJetsAndDPS DMuWJmix(lepSelection+"_13TeV_WJetsALL_MIX_UNFOLDING_dR_5311",    muLumi*36703.   *1000/76102995.,1, 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuWJmix.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //    }
        //    if ( (lepSelection.find("SE") == -1 && lepSelection.find("SMu") == -1 )  ) {
        //        //--- the following Bg are not used for WJets analysis ---//
        //        ZJetsAndDPS DMuTTrew(lepSelection+"_13TeV_TTJets_dR_5311_TopReweighting",             muLumi*245.           *1000/6923652., 1, 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuTTrew.Loop(1, 1, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //        ZJetsAndDPS DMuZZ(lepSelection+"_13TeV_ZZJets2L2Nu_dR_5311",        muLumi*17.654*0.04039 *1000/954911.,  1, 1, !doDataEff, bgSyst[i], bgDir[i], bgScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuZZ.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //        ZJetsAndDPS DMuWW(lepSelection+"_13TeV_WWJets2L2Nu_dR_5311",        muLumi*54.838*0.10608 *1000/1933235., 1, 1, !doDataEff, bgSyst[i], bgDir[i], bgScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuWW.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //        ZJetsAndDPS DMuZZ1(lepSelection+"_13TeV_ZZJets2L2Q_dR_5311",        muLumi*17.654*0.14118 *1000/1936727., 1, 1, !doDataEff, bgSyst[i], bgDir[i], bgScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuZZ1.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //        ZJetsAndDPS DMuZZ2(lepSelection+"_13TeV_ZZJets4L_dR_5311",          muLumi*17.654*0.010196*1000/4807893., 1, 1, !doDataEff, bgSyst[i], bgDir[i], bgScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuZZ2.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //        ZJetsAndDPS DMuWZ(lepSelection+"_13TeV_WZJets3LNu_dR_5311",         muLumi*33.21 *0.032887*1000/1995334., 1, 1, !doDataEff, bgSyst[i], bgDir[i], bgScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
        //        DMuWZ.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //        DMuWZ1.Loop(1, 0, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
        //    }
        }
    }
    
    // DYJets_MIX_UNFOLDING_dR_5311_Inf3
    if (doWhat == 3 || doWhat == 100 ){
        int doGen = 0 ;
        if ( lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0 ) doGen = 1;
        
        for (unsigned int i(0); i < NSystDYJets; i++){
            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && dySyst[i] == 4) continue; // jet smearing part -- not done for SMu ---
            if ((lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0) && dySyst[i] == 3) continue; // xsec -- not done for Z+jets ---
            if (dySyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuDYMix(lepSelection+"_13TeV_DYJets50toInf_dR_5311_List", muLumi * 5765.4 , 1., 1, !doDataEff, dySyst[i], dyDir[i], dyScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, 0);
            DMuDYMix.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);
        }
        
  //      if ( (lepSelection.find("SE") == -1 && lepSelection.find("SMu") == -1 )  ) {
  //          // these files are not needed for W+jets
  //          ZJetsAndDPS DMuDYTauS(lepSelection+"_13TeV_DYJets_FromTau_UNFOLDING_dR_5311_Inf3", muLumi*3531.8*1000/30459503., 1, 1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
  //          DMuDYTauS.Loop(1, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
  //          ZJetsAndDPS DMuDY(lepSelection+"_13TeV_DYJets_UNFOLDING_dR_5311_Inf3",  muLumi*3531.8*1000/30459503., 1., 1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
  //          DMuDY.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
  //          
  //          // scale files
  //          ZJetsAndDPS DMuDYscaleUp(lepSelection+"_13TeV_DYJets_UNFOLDING_dR_5311_Inf3_scaleUp",  muLumi*3531.8*1000/2170270., 1.,  1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
  //          DMuDYscaleUp.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
  //          ZJetsAndDPS DMuDYscaleDown(lepSelection+"_13TeV_DYJets_UNFOLDING_dR_5311_Inf3_scaleUp",  muLumi*3531.8*1000/1934901.,    1., 1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
  //          DMuDYscaleDown.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
  //      }
        
    }
        
    // W+jets inclusive signal samples 
//    if ( doWhat == 4 || doWhat == 100 ){
//        int doGen = 0 ;
//        if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
//        
//        for (unsigned int i(0); i < NSystWJets; i++){
//            if (!doGen ) continue;
//            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
//            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
//            
//            //ZJetsAndDPS DMuWJMix(lepSelection+"_13TeV_WJetsALL_MIX_UNFOLDING_dR_5311_List", muLumi* 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
//            //DMuWJMix.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
// 
//            ZJetsAndDPS DMuWJFxFx(lepSelection+"_13TeV_WJets_FxFx_dR_5311_List", muLumi* 60290.0 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
//            DMuWJFxFx.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
//
//            ZJetsAndDPS DMuWJMLM(lepSelection+"_13TeV_WJets_MLM_dR_5311_List", muLumi* 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
//            DMuWJMLM.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //open if you want to run on wjets by mg lo mlm sample as well
//
//            //ZJetsAndDPS DMuWJ(lepSelection+"_8TeV_WJetsALL_UNFOLDING_dR_5311",  muLumi*36703.       *1000/76102995., 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
//            //DMuWJ.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
//        }
//    }
 
    //W+jets inclusive NLO-FxFx sample   
    if (doWhat == 41 || doWhat == 100){
        int doGen = 0;
        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0) && lepSelection.find("SMuE") == -1)  doGen = 1 ;

        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

            ZJetsAndDPS DMuWJFxFx(lepSelection+"_13TeV_WJets_FxFx_dR_5311_List", muLumi* 60290.0 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, 0);
            DMuWJFxFx.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);
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

           ZJetsAndDPS DMuWJMLM(lepSelection+"_13TeV_WJets_MLM_dR_5311_List", muLumi* 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, 0);
           DMuWJMLM.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy);
       }
     }    
    
    // W+jets FxFx W pT-binned signal samples 
    if (doWhat == 5 || doWhat == 100 ){
        int doGen = 0 ;
        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuWJFxFx_Wpt1(lepSelection+"_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List", muLumi* 56306.4 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin, ZEtaMax, 0);
            DMuWJFxFx_Wpt1.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
             
            ZJetsAndDPS DMuWJFxFx_Wpt2(lepSelection+"_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List", muLumi* 3241.33 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_Wpt2.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
             
            ZJetsAndDPS DMuWJFxFx_Wpt3(lepSelection+"_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List", muLumi* 677.82 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_Wpt3.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
             
            ZJetsAndDPS DMuWJFxFx_Wpt4(lepSelection+"_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List", muLumi* 24.083 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_Wpt4.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
             
            ZJetsAndDPS DMuWJFxFx_Wpt5(lepSelection+"_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List", muLumi* 3.0563 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_Wpt5.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
             
            ZJetsAndDPS DMuWJFxFx_Wpt6(lepSelection+"_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List", muLumi* 0.4602 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_Wpt6.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
            
        }
    }
    
    // W+jets FxFx jet-binned signal samples 
    if ( doWhat == 6 || doWhat == 100 ){
        int doGen = 0 ;
        if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
        
        for (unsigned int i(0); i < NSystWJets; i++){
            if (!doGen ) continue;
            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
            
            ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", muLumi* 49264.92 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_jet0.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
             
            ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", muLumi* 8280.36 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_jet1.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
             
            ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", muLumi* 3118.08 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
            DMuWJFxFx_jet2.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); //FxFx with NLO normalization
         
        }
    }
    
 //   // Sherpa2
 //   if ( doWhat == 51){
 //       
 //       // this is setup for sherpa NLO
 //       ZJetsAndDPS DESherpaTest2NLO("SMu_13TeV_WToLNu_Sherpa2jNLO4jLO_v5",  muLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       DESherpaTest2NLO.Loop(0, 1, 0, 0, 0);
 //       
 //       // this is setup for sherpa NLO
 //       //ZJetsAndDPS DESherpaTest2NLO("DE_8TeV_DY_Sherpa_2NLO4_HepMC_dR_Full_List",  eleLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //DESherpaTest2NLO.Loop(0, 1, 0, 0, 0);
 //       
 //       
 //       // ZJetsAndDPS DESherpaTest1NLO("DE_8TeV_DY_Sherpa_1NLO4_HepMC_dR_Full_List",  eleLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //DESherpaTest1NLO.Loop(0, 1, 0, 0, 0);
 //       
 //       //ZJetsAndDPS DESherpaTest1NLOScaleDown("DE_8TeV_DY_Sherpa_1NLO4_scaleDown_HepMC_dR_Full_List",  eleLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //	DESherpaTest1NLOScaleDown.Loop(0, 1, 0, 0);
 //       
 //       //ZJetsAndDPS DESherpaTest1NLOScaleUp("DE_8TeV_DY_Sherpa_1NLO4_scaleUp_HepMC_dR_Full_List",  eleLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //      DESherpaTest1NLOScaleUp.Loop(0, 1, 0, 0);
 //       
 //       
 //       
 //       //ZJetsAndDPS DMuSherpa(lepSelection+"_DYJets_Sherpa_mcEveWeight",   muLumi * 3531.8           * 1000 / 30459503.,    1.,  0,  0,  0,  0,  1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //ZJetsAndDPS DESherpALL("DE_8TeV_Sherpa_HepMC_Z2jetNLO4jetLO_multithread_ALL_dR",  eleLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //DESherpALL.Loop(0, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
 //       
 //       
 //       //ZJetsAndDPS DMuPowMiNLO("DMu_8TeV_DYJets_PowhegZ2jMiNLO_dR_GEN_Cern",                             muLumi * 1.            * 1000 / 1964662.,    1.013,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax   );
 //       //DMuPowMiNLO.Loop(0, 1, 0, 0, 0);
 //       //	ZJetsAndDPS DEPow("DE_8TeV_DYJets_PowhegNLO1Jet_dR_GEN",                             muLumi * 334.            * 1000 / 2948078.,    1.013,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax   );
 //       //	DEPow.Loop(0, 1, 0, 0);
 //       //	ZJetsAndDPS DEPowSU("DE_8TeV_DYJets_PowhegNLO1Jet_dR_ScaleUp_GEN",                             muLumi * 318.4            * 1000 / 5446372.,    1.013,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax   );
 //       //	DEPowSU.Loop(0, 1, 0, 0);
 //       //	ZJetsAndDPS DEPowSD("DE_8TeV_DYJets_PowhegNLO1Jet_dR_ScaleDown_GEN",                             muLumi * 357.2            * 1000 / 5856584.,    1.013,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax   );
 //       //	DEPowSD.Loop(0, 1, 0, 0);
 //       
 //       //ZJetsAndDPS DESherpaTest("DE_8TeV_DY_Sherpa_HepMC_dR_Full_List",  eleLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //ZJetsAndDPS DESherpaTest("Sherpa\/DE_8TeV_Sherpa_HepMC_num*",  eleLumi         * 1000.          , 1.,    0,   0,     0,    0,     1.,  jetPtMin,  jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax );
 //       //DESherpaTest.Loop(0, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy ); 
 //       
 //   }
    
    
 //   // skip this part
 //   if ( doWhat == 2  ){
 //       for (unsigned int i(3); i < 5; i++){
 //           ZJetsAndDPS DMuDYTau(lepSelection+"_13TeV_DYJets_FromTau_UNFOLDING_dR_5311_Inf3", muLumi*3531.8*1000/3045950, 1., 1,  !doDataEff, tauSyst[i], tauDir[i], tauScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
 //           //	  DMuDYTau.Loop(1, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
 //       }
 //   }
 //   if (doWhat == 222){ // for individual production
 //       ZJetsAndDPS DMuDY("DMu_13TeV_DYJets_UNFOLDING_dR_5311_Inf3",  muLumi*3531.8*1000/30459503., 1., 1, !doDataEff, 0, 0, 1,          jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
 //       //       DMuDY.Loop(1, 1,  doQCD,  doSSign, doInvMassCut,  doBJets, doPUStudy );
 //       //         DMuDY.Loop(1, 0,  1,  1, 1, 22, doBJets, 9 );
 //       ZJetsAndDPS DEDY("DE_13TeV_DYJets_UNFOLDING_dR_5311_Inf3",  muLumi*3531.8*1000/30459503., 1., 1, !doDataEff, 0, 0, 1,          jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
 //       ///          DEDY.Loop(1, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
 //       ZJetsAndDPS DMuTTrew(lepSelection+"_13TeV_TTJets_dR_5311_TopReweighting",             muLumi*245.           *1000/6923652., 1, 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
 //       DMuTTrew.Loop(1, 1, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
 //       
 //   }
 //   if (doWhat == -3) { // this is a testing flag
 //       ZJetsAndDPS DMuDYMixPDF(lepSelection+"_13TeV_DYJets_MIX_UNFOLDING_dR_5311_Inf3", muLumi*3531.8*1000/30459503., 1., 1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, -30, 1);
 //       //DMuDYMixPDF.Loop(0, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, 0, 0, 1, 0, "CT10nlo.LHgrid", 0);
 //       //DMuDYMixPDF.Loop(1, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, 0, 0, 0, 0);
 //       
 //       //ZJetsAndDPS DMuDYscaleUp(lepSelection+"_8TeV_DYJets_UNFOLDING_dR_5311_Inf3_scaleUp",  muLumi*3531.8*1000/(2*2170270.), 1., 1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
 //       //DMuDYscaleUp.Loop(1, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
 //       //ZJetsAndDPS DMuDYscaleDown(lepSelection+"_8TeV_DYJets_UNFOLDING_dR_5311_Inf3_scaleDown",  muLumi*3531.8*1000/(2*1934901.), 1., 1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax);
 //       //DMuDYscaleDown.Loop(1, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy );
 //       
 //   }
 //   if (doWhat == 90) {
 //       for (int pdfMember(1); pdfMember <= 5; pdfMember++) {
 //           ZJetsAndDPS DMuDYMixPDF(lepSelection+"_13TeV_DYJets_MIX_UNFOLDING_dR_5311_Inf3", muLumi*3531.8*1000/30459503., 1., 1, !doDataEff, 0, 0, 1, jetPtMin, jetPtMax, ZEtaMin,    ZEtaMax, -30);
 //           DMuDYMixPDF.Loop(0, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, 0, 0, 1, 0, "NNPDF20_as_0118_100.LHgrid", pdfMember);
 //           //DMuDYMixPDF.Loop(0, 1,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, 0, 0, 1, 0);
 //       }
 //   }
 //   // now produce files for pulls
 //   if ( doWhat == 1001 ){
 //       int doGen = 1 ;
 //       doDataEff = 1 ;
 //       int NPulls = 25 ;
 //       for (int loopPull = 0 ; loopPull < NPulls ; loopPull++){
 //           if ( lepSelection.find("DMu") == 0 ||  lepSelection.find("DE") == 0 ) {
 //               ZJetsAndDPS DMuDYMix(lepSel+"_13TeV_DYJets_MIX_UNFOLDING_dR_5311_Inf3", muLumi*3531.8*1000/30459503., 1., 1, !doDataEff, 0, 0, 1,   jetPtMin,         jetPtMax, ZPtMin , ZEtaMin,    ZEtaMax, 0);
 //               DMuDYMix.Loop(1, 1,  0,  0, 0, 0, -10, 0, 0, 1, 0, "", 0, loopPull, NPulls  );

 //           }
 //       }


 //   }
        
        
    //andrew -- now that we have a script that recompiles before all of the runDYJet() batch submissions are done, we ask the runDYJets.cc code to load the libraries instead of recompiling      
//    //--- clean the *_cc.d and *_cc.so files ---
//    string cmd = "if ls *_cc.d &> .ls_tmp.list; then rm *_cc.d; fi";
//    system(cmd.c_str());
//    cmd = "if ls *_cc.so &> .ls_tmp.list; then rm *_cc.so; fi";
//    system(cmd.c_str());
//    cmd = "if ls " + srcdir + "*_cc.d &> .ls_tmp.list; then rm " + srcdir + "*_cc.d; fi";
//    system(cmd.c_str());
//    cmd = "if ls " + srcdir + "*_cc.so &> .ls_tmp.list; then rm " + srcdir + "*_cc.so; fi";
//    system(cmd.c_str());
//    system("rm .ls_tmp.list");

}

