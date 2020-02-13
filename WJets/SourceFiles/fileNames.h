#ifndef __fileNames_H___
#define __fileNames_H___

// Directory of input root files
const string FILESDIRECTORY("HistoFiles/");

// Basic information on samples in common struct
struct processInfoStruct{
    string filename;
    double NEvents, xsec, xsecFrac, xsecError;
    int color;
    string legend;
};

// the only three methods of processInfoStruct currently used are...
// filename, legend, color

const processInfoStruct ProcessInfo[] = {
    //--  Name  ---------------------------- #events ---- xsec - branch - xsec error (%) -- color for plot -- name on legend
    {"Data_dR_5311_List",                          1,          1.,      1,         1,             kBlack,      " Data"},   // 0
    {"WJetsALL_MIX_UNFOLDING_dR_5311_List",        76102995.,  36864.,  1,         0.03,          kOrange,       " W(#mu#nu)+jets"},
    {"ZZ_dR_5311_List",	                          9799908.,   17.654,  1,         0.04,          kOrange,     " ZZ"},
    {"WZ_dR_5311_List",	                          10000283.,  33.21,   1,         0.05,          kRed+1,      " WZ"},	
    {"WW_dR_5311_List",	                          10000431.,  54.838,  1,         0.05,          kViolet+5,   " WW"},
    {"VV_dR_5311_List",	                          10000431.,  54.838,  1,         0.05,          kRed+1,   " VV"}, //5

    {"ST_s_channel_dR_5311_List",	                  259961.,    3.79,    1,         0.10,          kMagenta,    " STs"},     // 6
    {"ST_tW_top_channel_dR_5311_List",	                  3758227.,   56.4,    1,         0.10,          kMagenta,    " STtWtop"},
    {"ST_tW_antitop_channel_dR_5311_List",	              497658.,    11.1,    1,         0.10,          kMagenta,    " STtWantitop"},
    {"TTJets_dR_5311_List",		                  6923652.,   234.,    1,         0.10,          kBlue,       " t#bar{t}"}, // 9
    {"DYJets50toInf_dR_5311_List",     30459503.,  3531.8,  1,         0.04,          kGreen-8,   " DY+jets"}, // 10 /// up to this line files are set for W+jet s and TTbar
    {"Tbar_s_channel_dR_5311_List",                139974.,    1.76,    1,         0.10,          kMagenta,    " Tbars"},
    {"Tbar_t_channel_dR_5311_List",                1903681.,   30.7,    1,         0.10,          kMagenta,    " Tbart"},
    {"Tbar_tW_channel_dR_5311_List",               493460.,    11.1,    1,         0.10,          kMagenta,    " TbartW"},
  
    {"DYJets10to50_dR_5311_List",	              11707222.,  860.5,   1,         0.04,          kAzure-4,    " DY"},
    
    {"ZZJets2L2Nu_dR_5311_List",		              954911.,    17.654,  0.04039,   0.04,          kSpring+5,   " ZZJets2L2Nu"},
    {"ZZJets4L_dR_5311_List",		              4807893.,   17.654,  0.010196,  0.04,          kOrange,     " ZZJets4L"},
    {"ZZJets2L2Q_dR_5311_List",		              1936727.,   17.654,  0.14118,   0.04,          kOrange-3,   " ZZJets2L2Q"},
    {"WZJets3LNu_dR_5311_List",		              1995334.,   33.21,   0.032887,  0.04,          kRed+1,      " WZJets3LNu"},
    {"WZJets2L2Q_dR_5311_List",		              3215990.,   33.21,   0.068258,  0.04,          kPink-6,     " WZJets2L2Q"},
    {"WWJets2L2Nu_dR_5311_List",		              1933235.,   54.838,  0.10608 ,  0.04,          kViolet+5,   " WWJets2L2Nu"},
    {"Top_dR_5311_List",		                      1.,         1,       1,         0.04,          kMagenta,    " Single top"},  // 21

    {"DYJets_FromTau_UNFOLDING_dR_5311_List_Inf3", 30459503.,  3531.8,  1,         0.033,         kAzure+4,    " DYtautau"},
    {"DYJets10toInf3_dR_5311_List",	              1.,         1,       1,         0.04,          kGreen-8,   " DY+jets"},          // 23
    {"DataQCD_dR_5311_List",	                      1.,         1,       1,         0.04,          kGreen+3,     " QCD multijet"},         // 24

    //{"WJetsALL_MIX_UNFOLDING_dR_5311_List",        76102995.,  36864.,  1,         0.03,          kOrange,       " W(#mu#nu)+jets"},       // 25
    //andrew -- 2016 update (from MIX to FxFx), these xsec numbers don't matter?
    //{"WJets_FxFx_dR_5311_List",       76102995.,  36864.,  1,         0.03,          kOrange,       " W(#mu#nu)+jets"},       // 25
    {"WJets_FxFx_012J_dR_5311_List",       76102995.,  36864.,  1,         0.03,          kOrange,       " W(#mu#nu)+jets"},       // 25
    // {"WJets_FxFx_Wpt_dR_5311_List",       76102995.,  36864.,  1,         0.03,          kOrange,       " W(#mu#nu)+jets"},       // 25

    {"WJetsALL_MIX_dR_5311_List",                  76102995.,  36864.,  1,         0.03,          kPink,       " W(#mu#nu)+jets"},  // relative weight for mixed DY and WJ files are set inthe code
    {"WJetsALL_dR_5311_List",                      76102995.,  36864.,  1,         0.03,          kPink,       " W(#mu#nu)+jets"},
    {"DYJets_UNFOLDING_dR_5311_List_Inf3",         30459503.,  3531.8,  1,         0.04,          kAzure+10,   " DY+jets"}, /// up to this line files are set for W+jet s and TTbar

    //

    {"ttV_dR_5311_List",       1.,  1.,  1,         1,          kViolet+6,       " t#bar{t}V"},       // 29
};


//////////////////////////////////////////////////////////////////////////////////////////

// first element must point to the data
// last element must point to the MC Signal
const int DATAFILENAME(0);

// andrew -- for Run II W+jets
// main event selection, signal & backgrounds
const int NFILESTTBARWJETS(8);
const int FilesTTbarWJets[NFILESTTBARWJETS] = {0, 29, 5, 24, 21, 10, 9, 25};

// QCD BG turned off here
// used for ttbar studies
const int NFILESTTBARWJETS_NOQCD(7);
const int FilesTTbarWJets_NoQCD[NFILESTTBARWJETS_NOQCD] = {0, 29, 5, 21, 10, 9, 25};

// data, QCD BG turned off here
// used for btag efficiencies studies
const int NFILESTTBARWJETS_NOQCD_NODATA(6);
const int FilesTTbarWJets_NoQCD_NoData[NFILESTTBARWJETS_NOQCD_NODATA] = {29, 5, 21, 10, 9, 25};

// some testing - 2 sept 2019 - andrew
const int NFILESTTBARWJETS_NOQCD_NOTTBAR(6);
const int FilesTTbarWJets_NoQCD_NoTTBar[NFILESTTBARWJETS_NOQCD_NOTTBAR] = {0, 29, 5, 21, 10, 25};

// looking at QCD contribution
const int NFILESTTBARWJETS_DATAQCDVV(2);
const int FilesTTbarWJets_DataQCDVV[NFILESTTBARWJETS_DATAQCDVV] = {0, 24};

//////////////////////////////////////////////////////////////////////////////////////////


#endif
