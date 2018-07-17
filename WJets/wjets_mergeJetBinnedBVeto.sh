#!/bin/bash

cd HistoFiles
echo $PWD

# Merge jet-binned W+Jets signal samples

hadd SMu_13TeV_WJets_FxFx_012J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_0J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_1J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_2J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto.root

hadd SMu_13TeV_WJets_FxFx_012J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_WJets_FxFx_0J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_WJets_FxFx_1J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root SMu_13TeV_WJets_FxFx_2J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD1.root

hadd SMu_13TeV_WJets_FxFx_012J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_WJets_FxFx_0J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_WJets_FxFx_1J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root SMu_13TeV_WJets_FxFx_2J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD2.root

hadd SMu_13TeV_WJets_FxFx_012J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_WJets_FxFx_0J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_WJets_FxFx_1J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root SMu_13TeV_WJets_FxFx_2J_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BVeto_QCD3.root
