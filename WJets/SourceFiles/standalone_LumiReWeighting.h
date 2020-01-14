#ifndef _standalone_LumiReWeighting_H_
#define _standalone_LumiReWeighting_H_

#include "TH1.h"
#include "TFile.h"
#include <string>
#include "TH1.h"
#include "TFile.h"
#include <string>
#include <vector>
#include <TROOT.h>
#include <string>
#include <iostream>

using namespace std;

class standalone_LumiReWeighting {
 public:

  standalone_LumiReWeighting(int year = 2016, int mode = 0); // 0: central, -1: down, +1: up
  virtual ~standalone_LumiReWeighting();
  double weight(int numVertices);
  void weightOOT_init();

 protected:
  TH1F*      weights_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Here we put the arrays of data points that represent the pileup profiles in data and in the corresponding MC !!!
// The script will then normalize each the pileup profiles from data and MC, and then divide data/MC to get the weights dist.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// For 2016 data, All Eras --------------------------
// entries go from 0 vertices to 99 vertices

// Central value -- 69.2 mb
double Data2016Golden[100] = {
  238795.98503718118, 837538.4873224057, 2308410.187221423, 3124677.103652871, 4474095.1162428055, 5987854.0581500735, 6997757.969821492, 12887181.93524674, 35249762.074759796, 78674741.20603172, 
  176891869.93958282, 359985138.61673516, 602598513.5015937, 876277115.5996239, 1174127468.327614, 1488583992.0162895, 1758758520.9094644, 1943266103.961885, 2048547245.779879, 2101089655.107086, 
  2132451623.7336464, 2148886888.1450434, 2128858410.827557, 2062575735.4727747, 1962841661.896556, 1841845968.3399818, 1704117192.4511967, 1554507841.7759407, 1399475755.7263808, 1243518058.849184, 
  1088804114.165804, 937283829.2614226, 792019217.527644, 656689817.1166494, 534437649.489473, 427098596.45254654, 335080475.2111975, 257704006.48866978, 193735840.89168292, 141820119.90494114, 
  100664563.8273283, 69009834.44930163, 45537907.489921615, 28846400.068187505, 17505824.605293203, 10162433.424883349, 5637701.435195699, 2987253.7100131637, 1511992.9952906747, 731842.6044675907, 
  339821.1895073065, 152545.14112224057, 67404.76060182758, 30489.671792789683, 15152.097732143027, 8975.906923604369, 6496.152446463112, 5434.803787897476, 4889.956802605675, 4521.715814087792, 
  4208.4641672412135, 3909.7626755255747, 3614.2740284133693, 3320.72222241939, 3031.095639861233, 2748.2366785905765, 2474.9765260805257, 2213.8172129637005, 1966.8152147158328, 1735.5463725797435, 
  1521.1090879440671, 1324.1490594809602, 1144.8977881328822, 983.2201586562815, 838.6675611109142, 710.5335856105921, 597.9097516646505, 499.7391508105774, 414.86631303079787, 342.08204436375576, 
  280.1624067395204, 227.90140298963124, 184.1372753332624, 147.77261422324227, 117.78870098216248, 93.25467141811235, 73.33219166097734, 57.27638782669201, 44.433775812633016, 34.23790566318003, 
  26.203375901826135, 19.918796002930055, 15.039187884158551, 11.27822683451309, 8.400634164145696, 6.214952158443832, 4.566859382647435, 3.333122552000094, 2.416230595998896, 1.7397169535292134
};

// Systematic variation up -- 69.2 mb * (1+0.046)
// double Data2016Golden_up[100] = {};

// Systematic variation down -- 69.2 mb * (1-0.046)
// double Data2016Golden_down[100] = {};

// For 2016 MC: RunIISummer16MiniAODv3,
// pileup profile used is PUMoriond17_94X_mcRun2_asymptotic_v3
double MC_2016_RunIISummer16MiniAODv3[100] = {
  7452.1054463968285, 9555.887563480723, 23994.37353194471, 35915.116404197935, 44993.63373003245, 58610.96581445624, 99209.54112376105, 289142.3628070835, 526762.5298897822, 991044.4263724477, 
  2049254.7734905777, 3765738.6546575567, 6000745.461260786, 8363033.34563334, 10961490.96843843, 13767999.999768965, 16433155.902336435, 18435715.229311496, 20021383.42681115, 21463741.23554644, 
  22438687.309315544, 22850753.330253053, 22649510.13672858, 21985329.67090352, 20953563.40930449, 19461716.97440592, 17844522.973224625, 16040988.015172094, 14358751.364395162, 12597805.070909241, 
  11160593.081541259, 9724331.565745817, 8536180.400275828, 7459025.248379151, 6568356.230628695, 5821915.581814764, 5241610.545917776, 4715749.623185973, 4308484.820753838, 3925525.141758298, 
  3658406.270056405, 3403187.65603275, 3133525.3178225937, 2859098.850717221, 2545661.618561562, 2232274.953187751, 1987193.9381233272, 1675201.0532435295, 1391479.6282801016, 1109000.5367780519, 
  874724.6703240119, 656638.7674421158, 486589.25768303016, 348883.6860503871, 238688.62306540442, 145008.91486740482, 97810.64803626314, 64468.52854972701, 43778.36168717409, 28381.369868770504, 
  16049.634656948854, 10384.265837211027, 7658.792990899299, 7735.921296880419, 5225.570925843588, 6334.048306853649, 6012.209316179389, 5645.351796612423, 5967.465495415458, 5906.8873644838895, 
  5716.397968938369, 5146.698510950091, 5032.696733463985, 5212.356626425358, 6106.833297547656, 0.0, 0.0, 0.0, 0.0, 0.0, 
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// For 2017 data, All Eras --------------------------
// entries go from 0 vertices to 99 vertices

// Central value -- 69.2 mb
double Data2017Golden[100] = {
  258795.9844797146, 1084345.5097850803, 2011926.6131890407, 3778605.4626604784, 4055892.8026737953, 5877873.501293559, 6451599.601498123, 6833441.043881902, 9252971.713964242, 21823224.25370885, 
  43703362.90829493, 82653160.28896391, 131503661.1067244, 188569646.50508666, 266485912.87301922, 375265980.7435169, 523518489.22086996, 697806142.5645981, 871137767.455353, 1031126320.1137837, 
  1171612126.84883, 1285226082.7737958, 1372214019.1310866, 1439766227.4483762, 1497729932.605791, 1552615701.8893952, 1599705291.6600506, 1629287408.8259134, 1634893120.9349217, 1614647758.1767573, 
  1571086422.050321, 1508717466.127388, 1430645853.5115113, 1339206158.2971907, 1238967061.1563163, 1136636309.5830903, 1038144895.8739655, 947025949.0522391, 865558426.9132662, 796565233.76876, 
  743697142.572942, 710312745.8505831, 698077284.1111457, 705979763.0060952, 729706766.1749616, 761454801.3020145, 790509814.8799732, 805071906.1184887, 795131551.268046, 755182057.9526083, 
  685874304.0128254, 594252890.2502905, 490929975.8825336, 387195052.2455524, 292350746.3858193, 212137371.85685125, 148608048.47548652, 100994565.04313445, 66913078.53334446, 43422072.11634715, 
  27718353.7513125, 17473770.49319081, 10918042.587712254, 6784988.928476276, 4208095.005512156, 2613520.1235203124, 1630715.0697904504, 1025139.0937299866, 650703.2326218836, 417552.67297088855, 
  270918.1281409907, 177583.1515461191, 117411.49562933933, 78141.25555506989, 52235.65418774079, 34999.910854366026, 23462.541588273867, 15711.096051388959, 10495.253692708962, 6986.605319786724, 
  4630.603109100947, 3053.3639547939233, 2001.7317581875102, 1303.9778409877424, 843.6246771808393, 541.805464441386, 345.2811501015348, 218.26118984618086, 136.80764652537633, 85.00542980968123, 
  52.34469822398521, 31.936666402091813, 19.302382481913746, 11.554773698521458, 6.849765308507507, 4.0206539272263875, 2.336548489106792, 1.3442083507406772, 0.7654801324872996, 0.4314625385234751
};

// Systematic variation up -- 69.2 mb * (1+0.046)
// double Data2017Golden_up[100] = {};

// Systematic variation down -- 69.2 mb * (1-0.046)
// double Data2017Golden_down[100] = {};

// For 2017 MC: RunIIFall17MiniAODv2,
// pileup profile used is PU2017_12Apr2018_94X_mc2017_realistic_v14
double MC_2017_RunIIFall17MiniAODv2[100] = {
  14436.651898762524, 2927.2473896108345, 4923.097882527312, 15966.803943331824, 24282.84766381715, 37322.40421753814, 50228.90407173136, 56016.87050118915, 150620.1838654302, 154478.8281517354, 
  290595.8317686392, 580859.021788459, 996395.0944136698, 1623757.4326870826, 2504393.198511597, 3625861.5904813656, 4893891.936980968, 6265906.094161519, 7577113.339657881, 8788794.173907476, 
  9759775.438711341, 10491520.757764285, 11003190.295798307, 11474410.597175889, 11937514.439882275, 12336285.368366987, 12758341.219269058, 13025252.958521755, 12926125.71737357, 13077411.18473664, 
  13018866.236944422, 13108812.565825192, 12990924.330043592, 13082666.924367987, 12968637.33287269, 12650831.406050624, 12092525.49483212, 11518252.779670287, 11100787.384901924, 10796553.241431689, 
  9985173.487711376, 9079988.76082399, 8223103.615865181, 7462019.294566364, 6762739.810198194, 6169306.9303043615, 5648988.706801035, 5475615.827316358, 5251614.873661364, 5280554.7058086535, 
  5236246.824865908, 5238974.487206227, 5253145.0257059345, 5250949.590163726, 5232055.538830783, 5278625.383665501, 5196928.570155453, 5186217.505843468, 4936802.722578673, 4569366.646832748, 
  4074395.7245894624, 3508971.279946224, 2882473.8102197414, 2294496.255006547, 1784423.3973668593, 1553769.6087354783, 1337286.3586038041, 964528.01487677, 677790.827394436, 485058.19812846807, 
  332774.8055189408, 225398.04900003425, 152017.27921047175, 104249.92408000404, 84158.36245131149, 64399.44257143836, 43043.84229723204, 29472.058945399993, 18361.824534831598, 11442.876159387808, 
  7118.5334247354385, 5122.68293181896, 3259.8891384302474, 2794.1906900830695, 2461.548941263656, 1064.4535962221216, 532.2267981110608, 332.641748819413, 199.5850492916478, 
  199.5850492916478, 133.0566995277652, 66.5283497638826, 0.0, 133.0566995277652, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// For 2018 data, All Eras --------------------------
// entries go from 0 vertices to 99 vertices

// // Central value -- 69.2 mb
// double Data2018Golden[100] = {
// };

// // Systematic variation up -- 69.2 mb * (1+0.046)
// // double Data2018Golden_up[100] = {};

// // Systematic variation down -- 69.2 mb * (1-0.046)
// // double Data2018Golden_down[100] = {};

// // For 2018 MC: RunIIAutumn18MiniAOD,
// // pileup profile used is 102X_upgrade2018_realistic_v15
// double MC_2018_RunIIAutumn18MiniAOD[100] = {
// };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

standalone_LumiReWeighting::standalone_LumiReWeighting(int year, int mode){
  
  std::vector<float> MC_distr;
  std::vector<float> Lumi_distr;
  MC_distr.clear();
  Lumi_distr.clear();

  if (year != 2016 && year != 2017 && year != 2018){
    std::cout << "Select a compatible year for PU reweighting." << std::endl;
    std::cout << "Setting year to 2016." << std::endl;
    year = 2016;
  }
  switch (mode){
    case 0:
      std::cout << "Using central value " << std::endl;
      break;
    case 1:
      std::cout << "Using +1 sigma 5% value " << std::endl;
      break;
    case -1:
      std::cout << "Using -1 sigma 5% value " << std::endl;
      break;
    default:
      std::cout << "Using central value " << std::endl;
      break;
  } // end of switch

  Int_t NBins = 100;
  if (year == 2016) NBins = 100;
  else if (year == 2017) NBins = 100;
  else NBins = 100;

  for(int i(0); i < NBins; ++i) {
    // ----------------------- 2016 -----------------------
    if (year == 2016){
	    switch (mode){
	      case 0:
          Lumi_distr.push_back(Data2016Golden[i]);
          break;
	      case 1:
          // Lumi_distr.push_back(Data2016Golden_up[i]); // no syst var yet!
          Lumi_distr.push_back(Data2016Golden[i]);
	        break;
	      case -1:
          // Lumi_distr.push_back(Data2016Golden_down[i]); // no syst var yet!
          Lumi_distr.push_back(Data2016Golden[i]);
	        break;
        default:
	        abort();
	        break;
	    } // end of switch
	    MC_distr.push_back(MC_2016_RunIISummer16MiniAODv3[i]);
    }
    // ----------------------- 2017 -----------------------
    else if (year == 2017){
	    switch (mode){
	      case 0:
          Lumi_distr.push_back(Data2017Golden[i]);
          break;
	      case 1:
          // Lumi_distr.push_back(Data2017Golden_up[i]);
          Lumi_distr.push_back(Data2017Golden[i]); // no syst var yet!
	        break;
	      case -1:
          // Lumi_distr.push_back(Data2017Golden_down[i]);
          Lumi_distr.push_back(Data2017Golden[i]); // no syst var yet!
	        break;
        default:
	        abort();
	        break;
	    } // end of switch
	    MC_distr.push_back(MC_2017_RunIIFall17MiniAODv2[i]);
    }
    else if (year == 2018){
      // USING 2017 FOR NOW!!!!!
      switch (mode){
	      case 0:
          Lumi_distr.push_back(Data2017Golden[i]);
          break;
	      case 1:
          // Lumi_distr.push_back(Data2017Golden_up[i]);
          Lumi_distr.push_back(Data2017Golden[i]); // no syst var yet!
	        break;
	      case -1:
          // Lumi_distr.push_back(Data2017Golden_down[i]);
          Lumi_distr.push_back(Data2017Golden[i]); // no syst var yet!
	        break;
        default:
	        abort();
	        break;
	    } // end of switch
	    MC_distr.push_back(MC_2017_RunIIFall17MiniAODv2[i]);
    }
    else{
      std::cout << "Select a compatible year for PU reweighting." << std::endl;
    }
  } // end of loop over bins

  // Check that vectors are of same length
  if( MC_distr.size() != Lumi_distr.size() ){   
    std::cout << "MC_distr.size() = " << MC_distr.size() << std::endl;
    std::cout << "Lumi_distr.size() = " << Lumi_distr.size() << std::endl;
    std::cerr << "ERROR: standalone_LumiReWeighting: input vectors have different sizes. Quitting... \n";
  }

  // Make histograms
  weights_  = new TH1F(Form("lumiNum_%d",mode), Form("lumiNum_%d",mode), NBins, 0., float(NBins));
  TH1F* den = new TH1F(Form("lumiDenom_%d",mode), Form("lumiDenom_%d",mode), NBins, 0., float(NBins));

  // ROOT histos and C++ vectors are indexed differently...
  for(int ibin = 1; ibin < NBins+1; ++ibin ) {
    // first for pileup profile from data
    weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
    // then MC pileup profile in denominator
    den->SetBinContent(ibin, MC_distr[ibin-1]);
  }

  // std::cout << " ----- Pileup Profile: Data Input ----- " << std::endl;
  // for(int ibin = 1; ibin < NBins+1; ++ibin){
  //   // std::cout << "   " << ibin << ": " << weights_->GetBinContent(ibin) << std::endl;
  //   std::cout << weights_->GetBinContent(ibin) << std::endl;
  // }
  // std::cout << " ----- Pileup Profile: MC Input ----- " << std::endl;
  // for(int ibin = 1; ibin<NBins+1; ++ibin){
  //   // std::cout << "   " << ibin << ": " << den->GetBinContent(ibin) << std::endl;
  //   std::cout << den->GetBinContent(ibin) << std::endl;
  // }

  // Check integrals, make sure things are normalized
  // We do this because we only want to re-weight for the shape difference in the pileup dist.'s
  float dataIntegral = weights_->Integral();
  weights_->Scale(1.0/dataIntegral);
  float mcIntegral = den->Integral();
  den->Scale(1.0/mcIntegral);

  // Now divide to get the distribution of weights
  weights_->Divide(weights_, den, 1, 1);

  // Printing out weights_ to troubleshoot
  // Number of vertices goes from 0 to 99
  // std::cout << "\nPrinting values for pileup weights vector: " << std::endl;
  // for(int ibin = 1; ibin < NBins+1; ++ibin) {
  //     //std::cout << "   # vtx = " << ibin-1 << ": " << weights_->GetBinContent(ibin) << std::endl;
  //     std::cout << weights_->GetBinContent(ibin) << std::endl;
  // }
  // std::cout << "\nweights_->Integral() = " << weights_->Integral() << std::endl;
 
  
  // Checking the normalization of the weights distribution (?)
//  double inte = 0;
//  for(int ibin = 1; ibin < NBins+1; ++ibin){
//    inte += weights_->GetBinContent(ibin) * MC_distr[ibin-1];
//  }
//
//  std::cout << "PU weight normalisation: " << inte << "\n";


  delete den; // free heap memory
}

standalone_LumiReWeighting::~standalone_LumiReWeighting()
{
}

double standalone_LumiReWeighting::weight(int numVertices) {
  //int binTemp = weights_->GetXaxis()->FindBin(numVertices);
  //return weights_->GetBinContent(binTemp);
  
  // histogram bin 1 should be indexing #vtx = 0
  return weights_->GetBinContent(numVertices+1);
}

#endif
