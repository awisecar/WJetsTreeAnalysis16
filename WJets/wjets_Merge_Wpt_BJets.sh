#!/bin/bash

cd HistoFiles
echo $PWD

#Just need doQCD=0 if doing ttbar study

hadd SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BJets.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BJets.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BJets.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BJets.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BJets.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BJets.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_0_JetPtMin_30_VarWidth_BJets.root

