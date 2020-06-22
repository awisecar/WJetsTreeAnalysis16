#!/bin/bash

cd DataW_txt_2016
echo $PWD

# data
mv Bonzais-SingleMuon-all-VJetPruner-SMu.txt SMu_13TeV_Data_dR_5311_List.txt 

# signal 
mv Bonzais-WToLNu_0J_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt                           SMu_13TeV_WJets_FxFx_0J_dR_5311_List.txt
mv Bonzais-WToLNu_1J_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt                           SMu_13TeV_WJets_FxFx_1J_dR_5311_List.txt
mv Bonzais-WToLNu_2J_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt                           SMu_13TeV_WJets_FxFx_2J_dR_5311_List.txt
mv Bonzais-WJetsToLNu_Wpt-0To50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt   SMu_13TeV_WJets_FxFx_Wpt-0To50_dR_5311_List.txt
mv Bonzais-WjetsToLNu_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt  SMu_13TeV_WJets_FxFx_Wpt-50To100_dR_5311_List.txt
mv Bonzais-WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt SMu_13TeV_WJets_FxFx_Wpt-100To250_dR_5311_List.txt
mv Bonzais-WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt SMu_13TeV_WJets_FxFx_Wpt-250To400_dR_5311_List.txt
mv Bonzais-WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt SMu_13TeV_WJets_FxFx_Wpt-400To600_dR_5311_List.txt
mv Bonzais-WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt SMu_13TeV_WJets_FxFx_Wpt-600ToInf_dR_5311_List.txt
mv Bonzais-WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-all-VJetPruner-SMu.txt              SMu_13TeV_WJets_MLM_dR_5311_List.txt

# background
mv Bonzais-DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-all-VJetPruner-SMu.txt                             SMu_13TeV_DYJets50toInf_dR_5311_List.txt 
mv Bonzais-ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1-all-VJetPruner-SMu.txt                    SMu_13TeV_ST_s_channel_dR_5311_List.txt
mv Bonzais-ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1-all-VJetPruner-SMu.txt SMu_13TeV_ST_t_antitop_channel_dR_5311_List.txt
mv Bonzais-ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1-all-VJetPruner-SMu.txt     SMu_13TeV_ST_t_top_channel_dR_5311_List.txt
mv Bonzais-ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1-all-VJetPruner-SMu.txt                  SMu_13TeV_ST_tW_antitop_channel_dR_5311_List.txt
mv Bonzais-ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1-all-VJetPruner-SMu.txt                      SMu_13TeV_ST_tW_top_channel_dR_5311_List.txt
mv Bonzais-TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-all-VJetPruner-SMu.txt                                              SMu_13TeV_TTJets_dR_5311_List.txt
mv Bonzais-ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8-all-VJetPruner-SMu.txt                                        SMu_13TeV_ttH_non_bb_channel_dR_5311_List.txt
mv Bonzais-ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8-all-VJetPruner-SMu.txt                                           SMu_13TeV_ttH_bb_channel_dR_5311_List.txt
mv Bonzais-TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8-all-VJetPruner-SMu.txt                        SMu_13TeV_ttW_LNu_channel_dR_5311_List.txt
mv Bonzais-TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8-all-VJetPruner-SMu.txt                         SMu_13TeV_ttW_QQ_channel_dR_5311_List.txt
mv Bonzais-TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8-all-VJetPruner-SMu.txt                                SMu_13TeV_ttZ_LLNuNu_channel_dR_5311_List.txt
mv Bonzais-TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8-all-VJetPruner-SMu.txt                                         SMu_13TeV_ttZ_QQ_channel_dR_5311_List.txt
mv Bonzais-WWTo2L2Nu_13TeV-powheg-all-VJetPruner-SMu.txt                                                              SMu_13TeV_WW_dR_5311_List.txt
mv Bonzais-WZ_TuneCUETP8M1_13TeV-pythia8-all-VJetPruner-SMu.txt                                                       SMu_13TeV_WZ_dR_5311_List.txt
mv Bonzais-ZZ_TuneCUETP8M1_13TeV-pythia8-all-VJetPruner-SMu.txt                                                       SMu_13TeV_ZZ_dR_5311_List.txt


