#!/bin/bash

#Script is made to cd into the HistoFiles directory and merge all of the indnividual data files (made from breaking up big data list into smaller lists in order to speed up batch jobs) into one big data file, as required by further analysis steps

cd HistoFiles
echo $PWD

hadd SMu_13TeV_Data_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_01_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_02_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_03_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_04_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_05_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_06_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_07_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_08_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_09_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_10_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root

hadd SMu_13TeV_Data_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_01_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_02_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_03_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_04_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_05_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_06_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_07_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_08_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_09_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_Data_10_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root

hadd SMu_13TeV_Data_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_01_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_02_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_03_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_04_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_05_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_06_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_07_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_08_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_09_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_Data_10_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root

hadd SMu_13TeV_Data_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_01_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_02_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_03_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_04_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_05_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_06_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_07_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_08_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_09_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_Data_10_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root


#Now do systematics ---

#JES:
hadd SMu_13TeV_Data_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_01_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_02_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_03_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_04_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_05_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_06_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_07_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_08_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_09_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_10_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Down_JetPtMin_30_VarWidth_BVeto.root

hadd SMu_13TeV_Data_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_01_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_02_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_03_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_04_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_05_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_06_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_07_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_08_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_09_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_Data_10_dR_5311_List_EffiCorr_0_TrigCorr_1_Syst_2_Up_JetPtMin_30_VarWidth_BVeto.root
