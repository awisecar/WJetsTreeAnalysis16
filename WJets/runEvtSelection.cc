void runEvtSelection(int doWhat = 0, int doQCD = 0, int doSysRunning = 0, int year = 2016)
// ^^^turn on above line if doing the job submission with wjets_jobsub_Condor.py
// and comment out the below lines defining doWhat, doSysRunning, doQCD, year
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

    // dit bonjour
    welcomeMessage();
        
    //------------------

    //  int doWhat       = 223;
    // int doWhat       = 63;
    //    int doWhat       = 14;
    //   int doWhat       = 41;
    //    int doWhat       = 100;
                              // 100 - all samples; 200 - all MC samples 
                              // 10, 11, ... - individual data samples, 1 - background , 2 - tau ?, 3 - DY, 
                              // 41 - W+jets incl. NLO-FxFx, 42 - W+jets incl. LO-MLM
                              // 5 - W+jets FxFx W pT-binned, 6 - W+jets FxFx jet-binned,
                              // 51 - MC gen, 90 - PDF Syst., 1001 - do pull DY samples


        //  int doQCD        = 0;
                             // 0-3 : 4 combination between isolation/anti-isolation and MT cuts for QCD BG estimation

        
        //  int doSysRunning = 0;
        // int doSysRunning = 2;
                             // 0 - no syst running, 100 - all systematic runnings,
                             // 1 - PU, 2 - JES, 3 - XSEC, 4 - JER, 5 - LepSF,
                             // 6 - BtagSF, 7 - MES, 8 - MER, 9 - WB, 10 - RESP,
                             // 11 - L1Prefire

        
    // int doBJets      = -1; //normal btag veto
    int doBJets      = 0; //no btag veto
    // int doBJets      = 2; //ttbar SFs
                            // 0 - no information on b-jets will be used ;
                            // 1, 2 .. require at least 1, 2, .. ; use 2 for ttbar systmatics;
                            // -1, -2, .. veto the event if you have 1 or more, 2 or more .. b-jets ;
                            // 101 - require exactly 1 b-jet

                            
    //    int year        = 2016;
    //    int year        = 2017;
    //   int year        = 2018;
                           // 2016, 2017, or 2018 data/MC

    //------------------

    // printout about b-tag veto info
    // bTagVetoMessage(doBJets);

    string lepSelection = "SMu";

    double muLumi(0.);
    if (year == 2016){

        // INT. LUMINOSITY FOR 2016, FOR (HLT_IsoMu24 || HLT_IsoTkMu24) triggers
        muLumi = 35916.982; // 2016 data, legacy re-reco, ReReco_07Aug2017_Collisions16 json

        // NOMINAL INT. LUMINOSITY FOR 2016 ---
        // (5746.010 + 2572.903 + 4242.292 + 4025.228 + 3104.509 + 7575.579 + 8650.628) pb^-1 = 35917.150 pb^-1
        // muLumi = 35917.150; // 2016 data, legacy re-reco, ReReco_07Aug2017_Collisions16 json

    }
    else if (year == 2017){

        // INT. LUMINOSITY FOR 2017, FOR (HLT_IsoMu27) trigger
        muLumi = 41524.945; // 2017 data, EOY2017ReReco_Collisions17 json

        // NOMINAL INT. LUMINOSITY FOR 2017 ---
        // (4793.961 + 9631.214 + 4247.682 + 9313.642 + 13538.716) pb^-1 = 41525.217 pb^-1
        // muLumi = 41525.217; // 2017 data, EOY2017ReReco_Collisions17 json

        // REDUCED INT. LUMINOSITY FOR HLT_Mu27 PRESCALED TRIGGER ---
        // ( (27.264+8.951) + (42.011+4.128) + 21.282 + 32.486 + 48.814 ) pb^-1 = 184.939 pb^-1
        // muLumi = 184.939;

    }
    else{

        // INT. LUMINOSITY FOR 2018, FOR (HLT_IsoMu24) trigger
        muLumi = 59717.255; // 2018 data, 17SeptEarlyReReco2018ABC_PromptEraD_Collisions18 json

        // NOMINAL INT. LUMINOSITY FOR 2018 ---
        // (14027.047 + 7060.622 + 6894.771 + 31742.979) pb^-1 = 59725.420 pb^-1
        // muLumi = 59725.420; // 2018 data, 17SeptEarlyReReco2018ABC_PromptEraD_Collisions18 json

    }

    std::cout << "\n Using " << muLumi << " pb^-1 as the integrated luminosity!\n" << std::endl;

    //switches
    bool doSSign  =  0;     // contribution of QCD to emu in TTbar
    bool doInvMassCut = 0;
    int  doPUStudy = -10 ;  // default int the ZJets
    bool doFlat   = 0;
    bool doRoch   = 0;
    bool doVarWidth = 1;
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

    const int NSystData(3), NSystMC(14); // set NSystMC to 5 if you want to do only PU, XSEC
    const int NSystWJets(18), NSystDYJets(16);

    // systematic variation tags for data --- 
    short dataSyst[NSystData]  = {0, 2, 2};
    short dataDir[NSystData]   = {0,-1, 1};

    //systematic variation tags for background (generic)
    short bgSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short bgDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float bgScale[NSystMC]     = {1, 1, 1,   0.,   0.,  1, 1, 1, 1, 1, 1, 1, 1, 1};
 
    // systematic variations for background (ttbar, diboson, single top, DY+jets) ---

    // --- DY+jets
    short dySyst[NSystDYJets]  = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 4, 4, 8, 11, 11};
    short dyDir[NSystDYJets]   = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1,-1, 1, 1, -1, 1};
    float dyScale[NSystDYJets] = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    // --- ttBar
    short ttSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short ttDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float ttScale[NSystMC]     = {1, 1, 1,  0.07,  0.07,  1, 1, 1, 1, 1, 1, 1, 1, 1};

    // --- single top
    short tcsSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short tcsDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float tcsScale[NSystMC]    = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    short tctSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short tctDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float tctScale[NSystMC]    = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    short tcwSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short tcwDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float tcwScale[NSystMC]    = {1, 1, 1, 0.07, 0.07,  1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    // --- diboson (WW, WZ, ZZ)
    short wwSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short wwDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float wwScale[NSystMC]     = {1, 1, 1, 0.03, 0.03,  1, 1, 1, 1, 1, 1, 1, 1, 1};

    short wzSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short wzDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float wzScale[NSystMC]     = {1, 1, 1, 0.04, 0.04,  1, 1, 1, 1, 1, 1, 1, 1, 1};

    short zzSyst[NSystMC]      = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short zzDir[NSystMC]       = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float zzScale[NSystMC]     = {1, 1, 1, 0.03, 0.03,  1, 1, 1, 1, 1, 1, 1, 1, 1};
        
    // --- ttV (ttW, ttZ, ttH)
    short ttWSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11};
    short ttWDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float ttWScale[NSystMC]    = {1, 1, 1, 0.01, 0.01,  1, 1, 1, 1, 1, 1, 1, 1, 1};

    short ttZSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11}; 
    short ttZDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1}; 
    float ttZScale[NSystMC]    = {1, 1, 1, 0.01, 0.01,  1, 1, 1, 1, 1, 1, 1, 1, 1};

    short ttHSyst[NSystMC]     = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 8, 11, 11}; 
    short ttHDir[NSystMC]      = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1, 1, -1, 1};
    float ttHScale[NSystMC]    = {1, 1, 1, 0.15, 0.15,  1, 1, 1, 1, 1, 1, 1, 1, 1};

    // systematic variations for signal (W+jets) ---
    short wjSyst[NSystWJets]   = {0, 1, 1,    3,    3,  5, 5, 6, 6, 7, 7, 4, 4, 8, 9, 10, 11, 11};
    short wjDir[NSystWJets]    = {0,-1, 1,   -1,    1, -1, 1,-1, 1,-1, 1,-1, 1, 1, 1,  1, -1, 1};
    float wjScale[NSystWJets]  = {1, 1, 1, 1., 1.,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  1, 1, 1}; // we do not systematically vary the BG on the signal
        
    //------------------------------------

    if (year == 2016){
       std::cout << " >>>>>>>>>> Doing 2016 data/MC! >>>>>>>>>>>>>>>>>" << std::endl;
       // Data
       if (doWhat == 10 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata0(lepSelection+"_13TeV_Data_dR_5311_List_0", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata0.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }

       if (doWhat == 11 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                   
               ZJetsAndDPS DMudata1(lepSelection+"_13TeV_Data_dR_5311_List_1", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata1.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }

       if (doWhat == 12 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                   
                ZJetsAndDPS DMudata2(lepSelection+"_13TeV_Data_dR_5311_List_2", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata2.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }

       if (doWhat == 13 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                   
               ZJetsAndDPS DMudata3(lepSelection+"_13TeV_Data_dR_5311_List_3", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata3.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }

       if (doWhat == 14 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMudata4(lepSelection+"_13TeV_Data_dR_5311_List_4", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata4.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);      
           }
       }

       if (doWhat == 15 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
                ZJetsAndDPS DMudata5(lepSelection+"_13TeV_Data_dR_5311_List_5", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata5.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }

       if (doWhat == 16 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMudata6(lepSelection+"_13TeV_Data_dR_5311_List_6", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata6.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);                   
           }
       }

       if (doWhat == 17 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMudata7(lepSelection+"_13TeV_Data_dR_5311_List_7", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata7.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }

       if (doWhat == 18 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMudata8(lepSelection+"_13TeV_Data_dR_5311_List_8", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata8.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }

       if (doWhat == 19 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100) continue;
                   
                ZJetsAndDPS DMudata9(lepSelection+"_13TeV_Data_dR_5311_List_9", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata9.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
       }
           
        // Background -- TTbar
        if (doWhat == 21 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuTT(lepSelection+"_13TeV_TTJets_dR_5311_List", year, muLumi * 831.76  , 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTT.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // Background -- WW
        if (doWhat == 221 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWW2L2Nu(lepSelection+"_13TeV_WWTo2L2Nu_dR_5311_List", year, muLumi * 12.178  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWW2L2Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWWLNuQQ(lepSelection+"_13TeV_WWToLNuQQ_dR_5311_List", year, muLumi * 49.997  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWWLNuQQ.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- WZ
        if (doWhat == 222 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWZ2L2Q(lepSelection+"_13TeV_WZTo2L2Q_dR_5311_List", year, muLumi * 5.595 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ2L2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWZ1L3Nu(lepSelection+"_13TeV_WZTo1L3Nu_dR_5311_List", year, muLumi * 3.033 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ1L3Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWZ1L1Nu2Q(lepSelection+"_13TeV_WZTo1L1Nu2Q_dR_5311_List", year, muLumi * 10.71 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ1L1Nu2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ZZ
        if (doWhat == 223 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuZZ2L2Q(lepSelection+"_13TeV_ZZTo2L2Q_dR_5311_List", year, muLumi * 3.22 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuZZ2L2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuZZ2L2Nu(lepSelection+"_13TeV_ZZTo2L2Nu_dR_5311_List", year, muLumi * 0.564 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuZZ2L2Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ST_s,ST_t
        if (doWhat == 23 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
                if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuT1(lepSelection+"_13TeV_ST_s_channel_dR_5311_List", year, muLumi * 3.475  ,  1., 1, !doDataEff, tcsSyst[i], tcsDir[i], tcsScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT3(lepSelection+"_13TeV_ST_t_antitop_channel_dR_5311_List", year, muLumi * 80.95   ,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT3.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT2(lepSelection+"_13TeV_ST_t_top_channel_dR_5311_List", year, muLumi * 136.02  ,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ST_tW
        if (doWhat == 24 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuT5(lepSelection+"_13TeV_ST_tW_antitop_channel_dR_5311_List", year, muLumi * 35.85 ,  1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT5.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT4(lepSelection+"_13TeV_ST_tW_top_channel_dR_5311_List", year, muLumi * 35.85  , 1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT4.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ttW
        if (doWhat == 25 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTW1(lepSelection+"_13TeV_ttW_LNu_channel_dR_5311_List", year, muLumi * 0.179,  1., 1, !doDataEff, ttWSyst[i], ttWDir[i], ttWScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTW1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTW2(lepSelection+"_13TeV_ttW_QQ_channel_dR_5311_List", year, muLumi * 0.371, 1., 1, !doDataEff, ttWSyst[i], ttWDir[i], ttWScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTW2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ttZ
        if (doWhat == 26 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTZ1(lepSelection+"_13TeV_ttZ_LLNuNu_channel_dR_5311_List", year, muLumi * 0.26,  1., 1, !doDataEff, ttZSyst[i], ttZDir[i], ttZScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTZ1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTZ2(lepSelection+"_13TeV_ttZ_QQ_channel_dR_5311_List", year, muLumi * 0.60, 1., 1, !doDataEff, ttZSyst[i], ttZDir[i], ttZScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTZ2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ttH
        if (doWhat == 27 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTH1(lepSelection+"_13TeV_ttH_bb_channel_dR_5311_List", year, muLumi * 0.291,  1., 1, !doDataEff, ttHSyst[i], ttHDir[i], ttHScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTH1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTH2(lepSelection+"_13TeV_ttH_non_bb_channel_dR_5311_List", year, muLumi * 0.213, 1., 1, !doDataEff, ttHSyst[i], ttHDir[i], ttHScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTH2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // DY+Jets
        if (doWhat == 30 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0 ) doGen = 1;
           
           for (unsigned int i(0); i < NSystDYJets; i++){
               if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && dySyst[i] == 4) continue; // jet smearing part -- not done for SMu ---
               if ((lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0) && dySyst[i] == 3) continue; // xsec -- not done for Z+jets ---
               if (dySyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuDYMix(lepSelection+"_13TeV_DYJets50toInf_dR_5311_List", year, muLumi * 6077.22 , 1., 1, !doDataEff, dySyst[i], dyDir[i], dyScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuDYMix.Loop(1, doGen, year,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
           
        }
    
        // W+jets inclusive NLO-FxFx sample   
        if (doWhat == 41 || doWhat == 100 || doWhat == 200){
            int doGen = 0;
            if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0) && lepSelection.find("SMuE") == -1)  doGen = 1;

            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx(lepSelection+"_13TeV_WJets_FxFx_dR_5311_List", year, muLumi * 1.01848535 * 60410.0 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }

        // W+jets inclusive LO-MLM sample
        if (doWhat == 42 || doWhat == 100 || doWhat == 200){
            int doGen = 0;
            if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;

            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJMLM(lepSelection+"_13TeV_WJets_MLM_dR_5311_List", year, muLumi * 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJMLM.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }    
       
        // W+jets FxFx W pT-binned signal sample - 0 to 50 pT
        if (doWhat == 51 || doWhat == 100 || doWhat == 200){
            int doGen = 0 ;
            if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
            
                ZJetsAndDPS DMuWJFxFx_Wpt1(lepSelection+"_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List", year, muLumi * 1.01848535 * 56306.4 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_Wpt1.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
            }
        }

        // W+jets FxFx W pT-binned signal sample - 50 to 100 pT
        if (doWhat == 52 || doWhat == 100 || doWhat == 200){
            int doGen = 0 ;
            if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx_Wpt2(lepSelection+"_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List", year, muLumi * 1.01848535 * 3241.33 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_Wpt2.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
            }
        }

        // W+jets FxFx W pT-binned signal sample - 100 to 250 pT
        if (doWhat == 53 || doWhat == 100 || doWhat == 200){
            int doGen = 0 ;
            if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx_Wpt3(lepSelection+"_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List", year, muLumi * 1.01848535 * 677.82 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_Wpt3.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
            }
        }

        // W+jets FxFx W pT-binned signal sample - 250 to Inf pT
        if (doWhat == 54 || doWhat == 100 || doWhat == 200){
            int doGen = 0 ;
            if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
        
            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx_Wpt4(lepSelection+"_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List", year, muLumi * 1.01848535 * 24.083 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_Wpt4.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

                ZJetsAndDPS DMuWJFxFx_Wpt5(lepSelection+"_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List", year, muLumi * 1.01848535 * 3.0563 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_Wpt5.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
            
                ZJetsAndDPS DMuWJFxFx_Wpt6(lepSelection+"_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List", year, muLumi * 1.01848535 * 0.4602 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_Wpt6.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
            }
        }
       
        // W+jets FxFx jet-binned signal sample - W+0J
        if (doWhat == 61 || doWhat == 100 || doWhat == 200){
            int doGen = 0 ;
            if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
            
            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
                ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", year, muLumi * 1.01848535 * 49264.92 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_jet0.Loop(1, doGen, year,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
            }
        }

        // W+jets FxFx jet-binned signal sample - W+1J
        if (doWhat == 62 || doWhat == 100 || doWhat == 200){
            int doGen = 0 ;
            if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
            
            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", year, muLumi * 1.01848535 * 8280.36 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_jet1.Loop(1, doGen, year,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

            }
        }

        // W+jets FxFx jet-binned signal sample - W+2J
        if (doWhat == 63 || doWhat == 100 || doWhat == 200){
            int doGen = 0 ;
            if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
            
            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", year, muLumi * 1.01848535 * 3118.08 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx_jet2.Loop(1, doGen, year,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
            }
        }

    }
    else if (year == 2017){
        std::cout << " >>>>>>>>>> Doing 2017 data/MC! >>>>>>>>>>>>>>>>>" << std::endl;

        // Data
        if (doWhat == 10 || doWhat == 100) {
            for (unsigned int i(0); i < NSystData; i++) {
                if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

                ZJetsAndDPS DMudata0(lepSelection+"_13TeV_Data_dR_5311_List_0", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata0.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }

        if (doWhat == 11 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata1(lepSelection+"_13TeV_Data_dR_5311_List_1", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata1.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 12 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata2(lepSelection+"_13TeV_Data_dR_5311_List_2", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata2.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 13 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata3(lepSelection+"_13TeV_Data_dR_5311_List_3", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata3.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 14 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata4(lepSelection+"_13TeV_Data_dR_5311_List_4", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata4.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 15 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata5(lepSelection+"_13TeV_Data_dR_5311_List_5", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata5.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 16 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata6(lepSelection+"_13TeV_Data_dR_5311_List_6", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata6.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 17 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata7(lepSelection+"_13TeV_Data_dR_5311_List_7", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata7.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 18 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata8(lepSelection+"_13TeV_Data_dR_5311_List_8", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata8.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 19 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata9(lepSelection+"_13TeV_Data_dR_5311_List_9", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata9.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // Background -- TTbar
        // /TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        if (doWhat == 211 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
               ZJetsAndDPS DMuTT1(lepSelection+"_13TeV_TT_FullHad_dR_5311_List", year, muLumi * 365.6, 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTT1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // /TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM
        if (doWhat == 212 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
               ZJetsAndDPS DMuTT2(lepSelection+"_13TeV_TT_SemiLep_dR_5311_List", year, muLumi * 371.65, 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTT2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // /TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        if (doWhat == 213 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
               ZJetsAndDPS DMuTT3(lepSelection+"_13TeV_TT_2L2Nu_dR_5311_List", year, muLumi * 94.45, 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTT3.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // Background -- WW
        if (doWhat == 221 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWW2L2Nu(lepSelection+"_13TeV_WWTo2L2Nu_dR_5311_List", year, muLumi * 12.178  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWW2L2Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWWLNuQQ(lepSelection+"_13TeV_WWToLNuQQ_dR_5311_List", year, muLumi * 49.997  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWWLNuQQ.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- WZ
        if (doWhat == 222 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWZ2L2Q(lepSelection+"_13TeV_WZTo2L2Q_dR_5311_List", year, muLumi * 5.595 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ2L2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWZ1L3Nu(lepSelection+"_13TeV_WZTo1L3Nu_dR_5311_List", year, muLumi * 3.033 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ1L3Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWZ1L1Nu2Q(lepSelection+"_13TeV_WZTo1L1Nu2Q_dR_5311_List", year, muLumi * 10.71 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ1L1Nu2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ZZ
        if (doWhat == 223 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuZZ2L2Q(lepSelection+"_13TeV_ZZTo2L2Q_dR_5311_List", year, muLumi * 3.22 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuZZ2L2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuZZ2L2Nu(lepSelection+"_13TeV_ZZTo2L2Nu_dR_5311_List", year, muLumi * 0.564 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuZZ2L2Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ST_s,ST_t
        // /ST_s-channel_4f_leptonDecays_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        // /ST_t-channel_top_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        // /ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        if (doWhat == 23 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuT1(lepSelection+"_13TeV_ST_s_channel_dR_5311_List", year, muLumi * 3.475,  1., 1, !doDataEff, tcsSyst[i], tcsDir[i], tcsScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT3(lepSelection+"_13TeV_ST_t_antitop_channel_dR_5311_List", year, muLumi * 80.95,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT3.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT2(lepSelection+"_13TeV_ST_t_top_channel_dR_5311_List", year, muLumi * 136.02,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
               
           }
        }

        // Background -- ST_tW
        // /ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        // /ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM
        if (doWhat == 24 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuT5(lepSelection+"_13TeV_ST_tW_antitop_channel_dR_5311_List", year, muLumi * 35.85,  1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT5.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT4(lepSelection+"_13TeV_ST_tW_top_channel_dR_5311_List", year, muLumi * 35.85, 1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT4.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ttW
        if (doWhat == 25 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTW1(lepSelection+"_13TeV_ttW_LNu_channel_dR_5311_List", year, muLumi * 0.179,  1., 1, !doDataEff, ttWSyst[i], ttWDir[i], ttWScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTW1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTW2(lepSelection+"_13TeV_ttW_QQ_channel_dR_5311_List", year, muLumi * 0.371, 1., 1, !doDataEff, ttWSyst[i], ttWDir[i], ttWScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTW2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ttZ
        if (doWhat == 26 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTZ1(lepSelection+"_13TeV_ttZ_LLNuNu_channel_dR_5311_List", year, muLumi * 0.26,  1., 1, !doDataEff, ttZSyst[i], ttZDir[i], ttZScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTZ1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTZ2(lepSelection+"_13TeV_ttZ_QQ_channel_dR_5311_List", year, muLumi * 0.60, 1., 1, !doDataEff, ttZSyst[i], ttZDir[i], ttZScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTZ2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ttH
        if (doWhat == 27 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTH1(lepSelection+"_13TeV_ttH_bb_channel_dR_5311_List", year, muLumi * 0.291,  1., 1, !doDataEff, ttHSyst[i], ttHDir[i], ttHScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTH1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTH2(lepSelection+"_13TeV_ttH_non_bb_channel_dR_5311_List", year, muLumi * 0.213, 1., 1, !doDataEff, ttHSyst[i], ttHDir[i], ttHScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTH2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // DY+Jets
        // /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM
        if (doWhat == 30 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0 ) doGen = 1;
           
           for (unsigned int i(0); i < NSystDYJets; i++){
               if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && dySyst[i] == 4) continue; // jet smearing part -- not done for SMu ---
               if ((lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0) && dySyst[i] == 3) continue; // xsec -- not done for Z+jets ---
               if (dySyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuDYMix(lepSelection+"_13TeV_DYJets50toInf_dR_5311_List", year, muLumi * 6077.22, 1., 1, !doDataEff, dySyst[i], dyDir[i], dyScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuDYMix.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
           
        }

        // W+jets inclusive NLO-FxFx sample   
        if (doWhat == 41 || doWhat == 100 || doWhat == 200){
            int doGen = 0;
            if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0) && lepSelection.find("SMuE") == -1)  doGen = 1;

            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx(lepSelection+"_13TeV_WJets_FxFx_dR_5311_List", year, muLumi * 0.91462316 * 67270.0 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }

        // W+jets inclusive LO-MLM sample
        // ...
        if (doWhat == 42 || doWhat == 100 || doWhat == 200){
            int doGen = 0;
            if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;

            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJMLM(lepSelection+"_13TeV_WJets_MLM_dR_5311_List", year, muLumi * 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJMLM.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }   
        
        // W+jets FxFx W pT-binned signal sample - 0 to 50 pT
        // ...
        if ( doWhat == 51 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt1(lepSelection+"_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List", year, muLumi * 0.91462316 * 62169.3, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt1.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx W pT-binned signal sample - 50 to 100 pT
        // /WJetsToLNu_Pt-50To100_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        if ( doWhat == 52 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt2(lepSelection+"_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List", year, muLumi * 0.91462316 * 3582.58, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt2.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx W pT-binned signal sample - 100 to 250 pT
        // /WJetsToLNu_Pt-100To250_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        if ( doWhat == 53 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt3(lepSelection+"_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List", year, muLumi * 0.91462316 * 770.99, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt3.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx W pT-binned signal sample - 250 to Inf pT
        // /WJetsToLNu_Pt-250To400_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        // /WJetsToLNu_Pt-400To600_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        // /WJetsToLNu_Pt-600ToInf_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        if ( doWhat == 54 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt4(lepSelection+"_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List", year, muLumi * 0.91462316 * 27.987, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt4.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

               ZJetsAndDPS DMuWJFxFx_Wpt5(lepSelection+"_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List", year, muLumi * 0.91462316 * 3.5841, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt5.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

               ZJetsAndDPS DMuWJFxFx_Wpt6(lepSelection+"_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List", year, muLumi * 0.91462316 * 0.5476, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt6.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

           }
        }

        // W+jets FxFx jet-binned signal sample - W+0J
        // /WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM
        if ( doWhat == 61 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", year, muLumi * 0.91462316 * 54549.40, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet0.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx jet-binned signal sample - W+1J
        // /WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM
        if ( doWhat == 62 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
            for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", year, muLumi * 0.91462316 * 8822.50, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet1.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

            }
        }

        // W+jets FxFx jet-binned signal sample - W+2J
        // /WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM
        if ( doWhat == 63 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", year, muLumi * 0.91462316 * 3312.48, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet2.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

    }
    else{
        std::cout << " >>>>>>>>> Doing 2018 data/MC! >>>>>>>>>>" << std::endl;

        // Data
        if (doWhat == 10 || doWhat == 100) {
            for (unsigned int i(0); i < NSystData; i++) {
                if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

                ZJetsAndDPS DMudata0(lepSelection+"_13TeV_Data_dR_5311_List_0", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMudata0.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }

        if (doWhat == 11 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata1(lepSelection+"_13TeV_Data_dR_5311_List_1", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata1.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 12 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata2(lepSelection+"_13TeV_Data_dR_5311_List_2", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata2.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 13 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata3(lepSelection+"_13TeV_Data_dR_5311_List_3", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata3.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 14 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata4(lepSelection+"_13TeV_Data_dR_5311_List_4", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata4.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 15 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata5(lepSelection+"_13TeV_Data_dR_5311_List_5", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata5.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 16 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata6(lepSelection+"_13TeV_Data_dR_5311_List_6", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata6.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 17 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata7(lepSelection+"_13TeV_Data_dR_5311_List_7", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata7.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 18 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata8(lepSelection+"_13TeV_Data_dR_5311_List_8", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata8.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        if (doWhat == 19 || doWhat == 100) {
           for (unsigned int i(0); i < NSystData; i++) {
               if (dataSyst[i] != doSysRunning && doSysRunning != 100)  continue;

               ZJetsAndDPS DMudata9(lepSelection+"_13TeV_Data_dR_5311_List_9", year, 1., 1., 1, doDataEff, dataSyst[i], dataDir[i], 1, jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMudata9.Loop(1, 0, year, doQCD, doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // Background -- TTbar
        // /TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM
        if (doWhat == 211 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
               ZJetsAndDPS DMuTT1(lepSelection+"_13TeV_TT_FullHad_dR_5311_List", year, muLumi * 365.6, 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTT1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // /TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext3-v2/MINIAODSIM
        if (doWhat == 212 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
               ZJetsAndDPS DMuTT2(lepSelection+"_13TeV_TT_SemiLep_dR_5311_List", year, muLumi * 371.65, 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTT2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // /TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        if (doWhat == 213 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;
               ZJetsAndDPS DMuTT3(lepSelection+"_13TeV_TT_2L2Nu_dR_5311_List", year, muLumi * 94.45, 1., 1, !doDataEff, ttSyst[i], ttDir[i], ttScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTT3.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
        }

        // Background -- WW
        if (doWhat == 221 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWW2L2Nu(lepSelection+"_13TeV_WWTo2L2Nu_dR_5311_List", year, muLumi * 12.178  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWW2L2Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWWLNuQQ(lepSelection+"_13TeV_WWToLNuQQ_dR_5311_List", year, muLumi * 49.997  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWWLNuQQ.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- WZ
        if (doWhat == 222 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWZ2L2Q(lepSelection+"_13TeV_WZTo2L2Q_dR_5311_List", year, muLumi * 5.595 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ2L2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWZ1L3Nu(lepSelection+"_13TeV_WZTo1L3Nu_dR_5311_List", year, muLumi * 3.033 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ1L3Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuWZ1L1Nu2Q(lepSelection+"_13TeV_WZTo1L1Nu2Q_dR_5311_List", year, muLumi * 10.71 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWZ1L1Nu2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ZZ
        if (doWhat == 223 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuZZ2L2Q(lepSelection+"_13TeV_ZZTo2L2Q_dR_5311_List", year, muLumi * 3.22 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuZZ2L2Q.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

                ZJetsAndDPS DMuZZ2L2Nu(lepSelection+"_13TeV_ZZTo2L2Nu_dR_5311_List", year, muLumi * 0.564 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuZZ2L2Nu.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ST_s,ST_t
        // /ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v4/MINIAODSIM
        // /ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        // /ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        if (doWhat == 23 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuT1(lepSelection+"_13TeV_ST_s_channel_dR_5311_List", year, muLumi * 3.475,  1., 1, !doDataEff, tcsSyst[i], tcsDir[i], tcsScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT3(lepSelection+"_13TeV_ST_t_antitop_channel_dR_5311_List", year, muLumi * 80.95,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT3.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT2(lepSelection+"_13TeV_ST_t_top_channel_dR_5311_List", year, muLumi * 136.02,  1., 1, !doDataEff, tctSyst[i], tctDir[i], tctScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
               
           }
        }

        // Background -- ST_tW
        // /ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM
        // /ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM
        if (doWhat == 24 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuT5(lepSelection+"_13TeV_ST_tW_antitop_channel_dR_5311_List", year, muLumi * 35.85,  1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT5.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuT4(lepSelection+"_13TeV_ST_tW_top_channel_dR_5311_List", year, muLumi * 35.85, 1., 1, !doDataEff, tcwSyst[i], tcwDir[i], tcwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuT4.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

           }
        }

        // Background -- ttW
        // /TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM
        // /TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        if (doWhat == 25 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTW1(lepSelection+"_13TeV_ttW_LNu_channel_dR_5311_List", year, muLumi * 0.179,  1., 1, !doDataEff, ttWSyst[i], ttWDir[i], ttWScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTW1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTW2(lepSelection+"_13TeV_ttW_QQ_channel_dR_5311_List", year, muLumi * 0.371, 1., 1, !doDataEff, ttWSyst[i], ttWDir[i], ttWScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTW2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ttZ
        // /TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM
        // /TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM
        if (doWhat == 26 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTZ1(lepSelection+"_13TeV_ttZ_LLNuNu_channel_dR_5311_List", year, muLumi * 0.26,  1., 1, !doDataEff, ttZSyst[i], ttZDir[i], ttZScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTZ1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTZ2(lepSelection+"_13TeV_ttZ_QQ_channel_dR_5311_List", year, muLumi * 0.60, 1., 1, !doDataEff, ttZSyst[i], ttZDir[i], ttZScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTZ2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // Background -- ttH
        // /ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM
        // /ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM
        if (doWhat == 27 || doWhat == 100 || doWhat == 200){
            for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuTTH1(lepSelection+"_13TeV_ttH_bb_channel_dR_5311_List", year, muLumi * 0.291,  1., 1, !doDataEff, ttHSyst[i], ttHDir[i], ttHScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTH1.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuTTH2(lepSelection+"_13TeV_ttH_non_bb_channel_dR_5311_List", year, muLumi * 0.213, 1., 1, !doDataEff, ttHSyst[i], ttHDir[i], ttHScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuTTH2.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

            }
        }

        // DY+Jets
        // /DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM
        if (doWhat == 30 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0 ) doGen = 1;
           
           for (unsigned int i(0); i < NSystDYJets; i++){
               if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && dySyst[i] == 4) continue; // jet smearing part -- not done for SMu ---
               if ((lepSelection.find("DMu") == 0 || lepSelection.find("DE") == 0) && dySyst[i] == 3) continue; // xsec -- not done for Z+jets ---
               if (dySyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuDYMix(lepSelection+"_13TeV_DYJets50toInf_dR_5311_List", year, muLumi * 6077.22, 1., 1, !doDataEff, dySyst[i], dyDir[i], dyScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuDYMix.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
           }
           
        }

        // W+jets inclusive NLO-FxFx sample   
        if (doWhat == 41 || doWhat == 100 || doWhat == 200){
            int doGen = 0;
            if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0) && lepSelection.find("SMuE") == -1)  doGen = 1;

            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJFxFx(lepSelection+"_13TeV_WJets_FxFx_dR_5311_List", year, muLumi * 0.92009421 * 66870.0 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJFxFx.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }

        // W+jets inclusive LO-MLM sample
        // ...
        if (doWhat == 42 || doWhat == 100 || doWhat == 200){
            int doGen = 0;
            if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;

            for (unsigned int i(0); i < NSystWJets; i++){
                if (!doGen ) continue;
                if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
                if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

                ZJetsAndDPS DMuWJMLM(lepSelection+"_13TeV_WJets_MLM_dR_5311_List", year, muLumi * 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
                DMuWJMLM.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
            }
        }   
        
        // W+jets FxFx W pT-binned signal sample - 0 to 50 pT
        // ...
        if ( doWhat == 51 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt1(lepSelection+"_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List", year, muLumi * 0.92009421 * 62478.6, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt1.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx W pT-binned signal sample - 50 to 100 pT
        // ...
        if ( doWhat == 52 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt2(lepSelection+"_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List", year, muLumi * 0.92009421 * 3574.23, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt2.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx W pT-binned signal sample - 100 to 250 pT
        // ...
        if ( doWhat == 53 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt3(lepSelection+"_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List", year, muLumi * 0.92009421 * 770.99, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt3.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx W pT-binned signal sample - 250 to Inf pT
        // ...
        // ...
        // ...
        if ( doWhat == 54 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_Wpt4(lepSelection+"_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List", year, muLumi * 0.92009421 * 28.052, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt4.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

               ZJetsAndDPS DMuWJFxFx_Wpt5(lepSelection+"_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List", year, muLumi * 0.92009421 * 3.5841, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt5.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

               ZJetsAndDPS DMuWJFxFx_Wpt6(lepSelection+"_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List", year, muLumi * 0.92009421 * 0.5490, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_Wpt6.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

           }
        }

        // W+jets FxFx jet-binned signal sample - W+0J
        // /WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        if ( doWhat == 61 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", year, muLumi * 0.92009421 * 54487.20, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet0.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }

        // W+jets FxFx jet-binned signal sample - W+1J
        // /WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        if ( doWhat == 62 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
            for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", year, muLumi * 0.92009421 * 9104.82, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet1.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

            }
        }

        // W+jets FxFx jet-binned signal sample - W+2J
        // /WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        if ( doWhat == 63 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", year, muLumi * 0.92009421 * 3376.80, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet2.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }
       
    }
    
}
