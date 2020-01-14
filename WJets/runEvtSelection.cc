void runEvtSelection(int doWhat = 0, int doQCD = 0, int doSysRunning = 0, int year = 2018)
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

    //  int doWhat       = 213;
    // int doWhat       = 100;
    //  int doWhat       = 63;
     //int doWhat       = 42;
                              // 100 - all samples; 200 - all MC samples 
                              // 10, 11, ... - individual data samples, 1 - background , 2 - tau ?, 3 - DY, 
                              // 41 - W+jets inc. NLO-FxFx, 42 - W+jets inc. LO-MLM
                              // 5 - W+jets FxFx W pT-binned, 6 - W+jets FxFx jet-binned,
                              // 51 - MC gen, 90 - PDF Syst., 1001 - do pull DY samples


    // int doQCD        = 0;
                             // 0-3 : 4 combination between isolation/anti-isolation and MT cuts for QCD BG estimation

        
    // int doSysRunning = 0;
                             // 0 - no syst running, 100 - all systematic runnings,
                             // 1 - PU, 2 - JES, 3 - XSEC, 4 - JER, 5 - LepSF,
                             // 6 - BtagSF, 7 - MES, 8 - MER, 9 - WB, 10 - RESP

        
    // int doBJets      = -1; //normal btag veto
    int doBJets      = 0; //no btag veto
    // int doBJets      = 2; //ttbar SFs
                            // 0 - no information on b-jets will be used ;
                            // 1, 2 .. require at least 1, 2, .. ; use 2 for ttbar systmatics;
                            // -1, -2, .. veto the event if you have 1 or more, 2 or more .. b-jets ;
                            // 101 - require exactly 1 b-jet

                            
    //   int year        = 2016;
    //  int year        = 2017;
    // int year        = 2018;
                           // 2016, 2017, or 2018 data/MC

    //------------------

    // printout about b-tag veto info
    bTagVetoMessage(doBJets);

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

    std::cout << "\n Using " << muLumi << " pb as the integrated luminosity!\n" << std::endl;

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

       // Background -- Diboson
       if (doWhat == 22 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWWInc(lepSelection+"_13TeV_WW_dR_5311_List", year, muLumi * 12.178  , 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWWInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuWZInc(lepSelection+"_13TeV_WZ_dR_5311_List", year, muLumi * 47.13 , 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWZInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuZZInc(lepSelection+"_13TeV_ZZ_dR_5311_List", year, muLumi * 16.523 ,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuZZInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

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
    
        //    //W+jets inclusive NLO-FxFx sample   
        //    if (doWhat == 41 || doWhat == 100 || doWhat == 200){
        //        int doGen = 0;
        //        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0) && lepSelection.find("SMuE") == -1)  doGen = 1 ;

        //        for (unsigned int i(0); i < NSystWJets; i++){
        //            if (!doGen ) continue;
        //            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
        //            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

        //            ZJetsAndDPS DMuWJFxFx(lepSelection+"_13TeV_WJets_FxFx_dR_5311_List", muLumi* 60290.0 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJFxFx.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
        //        }
        //    }

        //    //W+jets inclusive LO-MLM sample
        //    if (doWhat == 42 || doWhat == 100 || doWhat == 200){
        //    int doGen = 0;
        //    if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;

        //        for (unsigned int i(0); i < NSystWJets; i++){
        //            if (!doGen ) continue;
        //            if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue;
        //            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

        //            ZJetsAndDPS DMuWJMLM(lepSelection+"_13TeV_WJets_MLM_dR_5311_List", muLumi* 61526.7 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJMLM.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);
        //        }
        //    }    
       
        //    // W+jets FxFx W pT-binned signal sample - 0 to 50 pT
        //    if (doWhat == 51 || doWhat == 100 || doWhat == 200){
        //        int doGen = 0 ;
        //        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
            
        //        for (unsigned int i(0); i < NSystWJets; i++){
        //            if (!doGen ) continue;
        //            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
        //            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
                
        //            ZJetsAndDPS DMuWJFxFx_Wpt1(lepSelection+"_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List", muLumi* 56306.4 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJFxFx_Wpt1.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
        //        }
        //    }

        //    // W+jets FxFx W pT-binned signal sample - 50 to 100 pT
        //    if (doWhat == 52 || doWhat == 100 || doWhat == 200){
        //        int doGen = 0 ;
        //        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
            
        //        for (unsigned int i(0); i < NSystWJets; i++){
        //            if (!doGen ) continue;
        //            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
        //            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

        //            ZJetsAndDPS DMuWJFxFx_Wpt2(lepSelection+"_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List", muLumi* 3241.33 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJFxFx_Wpt2.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
        //        }
        //    }

        //    // W+jets FxFx W pT-binned signal sample - 100 to 250 pT
        //    if (doWhat == 53 || doWhat == 100 || doWhat == 200){
        //        int doGen = 0 ;
        //        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
            
        //        for (unsigned int i(0); i < NSystWJets; i++){
        //            if (!doGen ) continue;
        //            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
        //            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

        //            ZJetsAndDPS DMuWJFxFx_Wpt3(lepSelection+"_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List", muLumi* 677.82 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJFxFx_Wpt3.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
        //        }
        //    }

        //    // W+jets FxFx W pT-binned signal sample - 250 to Inf pT
        //    if (doWhat == 54 || doWhat == 100 || doWhat == 200){
        //        int doGen = 0 ;
        //        if ((lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1)  doGen = 1 ;
            
        //        for (unsigned int i(0); i < NSystWJets; i++){
        //            if (!doGen ) continue;
        //            if ((lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
        //            if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

        //            ZJetsAndDPS DMuWJFxFx_Wpt4(lepSelection+"_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List", muLumi* 24.083 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJFxFx_Wpt4.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

        //            ZJetsAndDPS DMuWJFxFx_Wpt5(lepSelection+"_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List", muLumi* 3.0563 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJFxFx_Wpt5.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
                
        //            ZJetsAndDPS DMuWJFxFx_Wpt6(lepSelection+"_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List", muLumi* 0.4602 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
        //            DMuWJFxFx_Wpt6.Loop(1, doGen,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
        //        }
        //    }

       
       // W+jets FxFx jet-binned signal sample - W+0J
       if ( doWhat == 61 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", year, muLumi* 50131.98 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet0.Loop(1, doGen, year,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
       }

       // W+jets FxFx jet-binned signal sample - W+1J
       if ( doWhat == 62 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", year, muLumi* 8426.09 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet1.Loop(1, doGen, year,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization

           }
       }

       // W+jets FxFx jet-binned signal sample - W+2J
       if ( doWhat == 63 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", year, muLumi* 3172.96 , 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet2.Loop(1, doGen, year,  doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
       }
    }

    else if (year == 2017){
        std::cout << " >>>>>>>>> Doing 2017 data/MC! >>>>>>>>>>" << std::endl;

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

       // Background -- Diboson
       // /WWTo2L2Nu_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM
       // /WZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
       // /ZZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM
        if (doWhat == 22 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWWInc(lepSelection+"_13TeV_WW_dR_5311_List", year, muLumi * 12.178, 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWWInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuWZInc(lepSelection+"_13TeV_WZ_dR_5311_List", year, muLumi * 47.13, 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWZInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuZZInc(lepSelection+"_13TeV_ZZ_dR_5311_List", year, muLumi * 16.523,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuZZInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

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


       // W+jets FxFx jet-binned signal sample - W+0J
       // /WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM
        if ( doWhat == 61 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", year, muLumi * 50131.98, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
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

               ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", year, muLumi * 8426.09, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
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

               ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", year, muLumi * 3172.96, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
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

        // Background -- Diboson
        // /WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        // /WZ_TuneCP5_PSweights_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        // /ZZ_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM
        if (doWhat == 22 || doWhat == 100 || doWhat == 200){
           for (unsigned int i(0); i < NSystMC; i++){
               if (bgSyst[i] != doSysRunning && doSysRunning != 100) continue;

               ZJetsAndDPS DMuWWInc(lepSelection+"_13TeV_WW_dR_5311_List", year, muLumi * 12.178, 1., 1, !doDataEff, wwSyst[i], wwDir[i], wwScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWWInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuWZInc(lepSelection+"_13TeV_WZ_dR_5311_List", year, muLumi * 47.13, 1., 1, !doDataEff, wzSyst[i], wzDir[i], wzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWZInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

               ZJetsAndDPS DMuZZInc(lepSelection+"_13TeV_ZZ_dR_5311_List", year, muLumi * 16.523,  1., 1, !doDataEff, zzSyst[i], zzDir[i], zzScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuZZInc.Loop(1, 0, year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth);

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

        // W+jets FxFx jet-binned signal sample - W+0J
        // /WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM
        if ( doWhat == 61 || doWhat == 100 || doWhat == 200){
           int doGen = 0 ;
           if ( (lepSelection.find("SE") == 0 || lepSelection.find("SMu") == 0 ) && lepSelection.find("SMuE") == -1 )  doGen = 1 ;
           
           for (unsigned int i(0); i < NSystWJets; i++){
               if (!doGen ) continue;
               if ( ( lepSelection.find("SMu") == 0 || lepSelection.find("SE") == 0 ) && wjSyst[i] == 3) continue; // xsec -- not done for SMu ---
               if (wjSyst[i] != doSysRunning && doSysRunning != 100) continue;
               
               ZJetsAndDPS DMuWJFxFx_jet0(lepSelection+"_13TeV_WJets_FxFx_0J_dR_5311_List", year, muLumi * 50131.98, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
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

               ZJetsAndDPS DMuWJFxFx_jet1(lepSelection+"_13TeV_WJets_FxFx_1J_dR_5311_List", year, muLumi * 8426.09, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
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

               ZJetsAndDPS DMuWJFxFx_jet2(lepSelection+"_13TeV_WJets_FxFx_2J_dR_5311_List", year, muLumi * 3172.96, 1., 1, !doDataEff, wjSyst[i], wjDir[i], wjScale[i], jetPtMin, jetPtMax, ZPtMin, ZEtaMin, ZEtaMax, METcut, jetEtaMin, jetEtaMax);
               DMuWJFxFx_jet2.Loop(1, doGen,  year, doQCD,  doSSign, doInvMassCut, doBJets, doPUStudy, doFlat, doRoch, doVarWidth); //FxFx with NLO normalization
           }
        }
       
    }

}
