#!/bin/bash

cd HistoFiles
echo $PWD

## Syst for WJets Signal is 1=PU, 4=JER, 5=LepSF, 6=BtagSF

## 1=PU
hadd -f SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Up_JetPtMin_30_VarWidth_BVeto.root

hadd -f SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_1_Down_JetPtMin_30_VarWidth_BVeto.root

## 4=JER
hadd -f SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Up_JetPtMin_30_VarWidth_BVeto.root

hadd -f SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_4_Down_JetPtMin_30_VarWidth_BVeto.root

## 5=LepSF
##hadd SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Up_JetPtMin_30_VarWidth_BVeto.root
##
##hadd SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_5_Down_JetPtMin_30_VarWidth_BVeto.root
##
#### 6=BtagSF
hadd -f SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Up_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Up_JetPtMin_30_VarWidth_BVeto.root

hadd -f SMu_13TeV_WJets_FxFx_Wpt_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Down_JetPtMin_30_VarWidth_BVeto.root SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List_EffiCorr_1_TrigCorr_1_Syst_6_Down_JetPtMin_30_VarWidth_BVeto.root
