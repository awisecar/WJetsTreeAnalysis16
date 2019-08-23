/**
   \class    standalone_LumiReWeighting standalone_LumiReWeighting.h "PhysicsTools/Utilities/interface/standalone_LumiReWeighting.h"
   \brief    Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data
   This class will trivially take two histograms:
   1. The generated "flat-to-N" distributions from a given processing (or any other generated input)
   2. A histogram generated from the "estimatePileup" macro here:
   https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc#How_to_use_script_estimatePileup
   and produce weights to convert the input distribution (1) to the latter (2).
   \author Salvatore Rappoccio, modified by Mike Hildreth
  
*/
#ifndef standalone_LumiReWeighting_h
#define standalone_LumiReWeighting_h
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

  //use a "year code": 2016, 2017, etc.
  standalone_LumiReWeighting(int year=2016,int mode=0); // 0: central, -1: down, +1: up
  virtual ~standalone_LumiReWeighting();
  double weight(int npv);
  void weightOOT_init();

 protected:

  TH1F*      weights_;
};

double Data2016Golden[125] = {238797, 837542.9, 2308427, 3124754, 4476191, 5995911, 7000896, 1.289165e+07, 3.526173e+07, 7.870123e+07, 1.769458e+08, 3.600895e+08, 6.027665e+08, 8.765194e+08, 1.174474e+09, 1.489059e+09, 1.759352e+09, 1.943926e+09, 2.049172e+09, 2.101582e+09, 2.132787e+09, 2.149099e+09, 2.128986e+09, 2.062649e+09, 1.962884e+09, 1.841872e+09, 1.704136e+09, 1.554523e+09, 1.399489e+09, 1.243533e+09, 1.088821e+09, 9.373048e+08, 7.920441e+08, 6.567177e+08, 5.344668e+08, 4.271268e+08, 3.351056e+08, 2.577246e+08, 1.937514e+08, 1.418309e+08, 1.006714e+08, 6.901386e+07, 4.554008e+07, 2.884748e+07, 1.750632e+07, 1.016264e+07, 5637781, 2987282, 1512002, 731845.4, 339822, 152545.4, 67404.82, 30489.69, 15152.11, 8975.911, 6496.155, 5434.805, 4889.958, 4521.716, 4208.464, 3909.763, 3614.274, 3320.722, 3031.096, 2748.237, 2474.977, 2213.817, 1966.815, 1735.546, 1521.109, 1324.149, 1144.898, 983.2202, 838.6676, 710.5336, 597.9098, 499.7392, 414.8663, 342.082, 280.1624, 227.9014, 184.1373, 147.7726, 117.7887, 93.25467, 73.33219, 57.27639, 44.43378, 34.23791, 26.20338, 19.9188, 15.03919, 11.27823, 8.400634, 6.214952, 4.566859, 3.333123, 2.416231, 1.739717, 1.244147, 0.8837248, 0.6234687, 0.4368819, 0.304064, 0.2101925, 0.1443179, 0.09841784, 0.06666193, 0.04484685, 0.02996649, 0.01988791, 0.01310968, 0.008583092, 0.005581416, 0.003604903, 0.002312551, 0.001473456, 0.0009324613, 0.0005861009, 0.0003658993, 0.000226881, 0.0001397275, 8.546985e-05, 5.192671e-05};
//with 69.2mb as x-sec (Metin, from pileup twiki)
double Data2016Golden_up[125] = {232683.1, 659469.5, 2183870, 2745817, 4071554, 5399909, 6385964, 9041697, 2.378256e+07, 5.400185e+07, 1.160124e+08, 2.460465e+08, 4.433764e+08, 6.801008e+08, 9.378589e+08, 1.218809e+09, 1.50169e+09, 1.730463e+09, 1.882647e+09, 1.968689e+09, 2.012668e+09, 2.040373e+09, 2.054739e+09, 2.036914e+09, 1.978317e+09, 1.889437e+09, 1.781252e+09, 1.657948e+09, 1.523143e+09, 1.382248e+09, 1.239674e+09, 1.097676e+09, 9.577106e+08, 8.219349e+08, 6.933169e+08, 5.748257e+08, 4.68639e+08, 3.757723e+08, 2.961771e+08, 2.290778e+08, 1.733823e+08, 1.279647e+08, 9.17574e+07, 6.370686e+07, 4.270414e+07, 2.757313e+07, 1.711861e+07, 1.020579e+07, 5837524, 3201829, 1684026, 849969, 412581.5, 193615.1, 88891.36, 40974.85, 19936.6, 11013.4, 7298.96, 5723.065, 4985.533, 4560.541, 4244.792, 3963.568, 3691.523, 3421.666, 3153.597, 2889.128, 2630.654, 2380.541, 2140.893, 1913.458, 1699.602, 1500.307, 1316.191, 1147.53, 994.2979, 856.2033, 732.7331, 623.1949, 526.7585, 442.4958, 369.417, 306.5032, 252.734, 207.1108, 168.6757, 136.5252, 109.8205, 87.7939, 69.75179, 55.07523, 43.21823, 33.70443, 26.1226, 20.12124, 15.40287, 11.71809, 8.85973, 6.657203, 4.971307, 3.689408, 2.721132, 1.99457, 1.452965, 1.051882, 0.7568067, 0.5411383, 0.3845357, 0.2715627, 0.1905936, 0.1329385, 0.09215045, 0.06348164, 0.04346127, 0.0295706, 0.01999499, 0.01343647, 0.008973298, 0.005955545, 0.003928198, 0.002574942, 0.001677427, 0.001085979, 0.0006987167};
//with 69.2mb + %4.6 as x-sec (Metin, from pileup twiki)
double Data2016Golden_dn[125] = {247408.7, 1069214, 2428245, 3566876, 4991685, 6592990, 8096912, 1.996092e+07, 5.191008e+07, 1.197673e+08, 2.727459e+08, 5.131262e+08, 8.022785e+08, 1.118695e+09, 1.462916e+09, 1.780146e+09, 2.005132e+09, 2.135028e+09, 2.198451e+09, 2.233931e+09, 2.252533e+09, 2.2297e+09, 2.154207e+09, 2.041692e+09, 1.905627e+09, 1.751113e+09, 1.584631e+09, 1.413642e+09, 1.242572e+09, 1.073682e+09, 9.098191e+08, 7.551503e+08, 6.139337e+08, 4.891657e+08, 3.820187e+08, 2.920884e+08, 2.180386e+08, 1.582419e+08, 1.111195e+08, 7.515025e+07, 4.87546e+07, 3.024635e+07, 1.790095e+07, 1.008984e+07, 5410208, 2758284, 1337443, 617797.6, 273138.8, 116942.4, 49859.98, 22465.64, 11743.57, 7638.636, 6018.585, 5280.96, 4837.555, 4483.996, 4154.032, 3829.045, 3506.22, 3187.679, 2876.816, 2577.045, 2291.368, 2022.216, 1771.414, 1540.186, 1329.193, 1138.586, 968.0713, 816.9839, 684.3606, 569.0134, 469.5982, 384.6778, 312.777, 252.4292, 202.2142, 160.7869, 126.8986, 99.40994, 77.29821, 59.65893, 45.70331, 34.75246, 26.22944, 19.64978, 14.61137, 10.78422, 7.900425, 5.744812, 4.146335, 2.970406, 2.112174, 1.490752, 1.044343, 0.7261761, 0.5011893, 0.3433382, 0.2334548, 0.157559, 0.1055464, 0.07017827, 0.04631491, 0.03033875, 0.01972572, 0.01272993, 0.008154115, 0.005184242, 0.00327153, 0.002049149, 0.001273951, 0.0007861199, 0.000481483, 0.0002927044, 0.0001766174, 0.0001057775, 6.287945e-05, 3.710048e-05, 2.172728e-05, 1.262951e-05, 7.286554e-06, 4.172663e-06, 2.371682e-06};
//with 69.2mb - %4.6 as x-sec (Metin, from pileup twiki)
double MC_2016_25ns_Moriond17[125] = {1.78653e-05, 2.56602e-05, 5.27857e-05, 8.88954e-05, 0.000109362, 0.000140973, 0.000240998, 0.00071209, 0.00130121, 0.00245255, 0.00502589, 0.00919534, 0.0146697, 0.0204126, 0.0267586, 0.0337697, 0.0401478, 0.0450159, 0.0490577, 0.0524855, 0.0548159, 0.0559937, 0.0554468, 0.0537687, 0.0512055, 0.0476713, 0.0435312, 0.0393107, 0.0349812, 0.0307413, 0.0272425, 0.0237115, 0.0208329, 0.0182459, 0.0160712, 0.0142498, 0.012804, 0.011571, 0.010547, 0.00959489, 0.00891718, 0.00829292, 0.0076195, 0.0069806, 0.0062025, 0.00546581, 0.00484127, 0.00407168, 0.00337681, 0.00269893, 0.00212473, 0.00160208, 0.00117884, 0.000859662, 0.000569085, 0.000365431, 0.000243565, 0.00015688, 9.88128e-05, 6.53783e-05, 3.73924e-05, 2.61382e-05, 2.0307e-05, 1.73032e-05, 1.435e-05, 1.36486e-05, 1.35555e-05, 1.37491e-05, 1.34255e-05, 1.33987e-05, 1.34061e-05, 1.34211e-05, 1.34177e-05, 1.32959e-05, 1.33287e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 

standalone_LumiReWeighting::standalone_LumiReWeighting(int year,int mode) {
  
  std::vector<float> MC_distr;
  std::vector<float> Lumi_distr;
  MC_distr.clear();
  Lumi_distr.clear();

  if (year != 2016 && year != 2017){
      std::cout << "Select a compatible year for PU reweighting." << std::endl;
      std::cout << "Setting year to 2016." << std::endl;
      year = 2016;
}
  switch (mode)
    {
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

  Int_t NBins = 60;
  if(year == 2016) NBins = 125;
  else if (year == 2017) NBins = 125;

  for(int i=0; i< NBins; ++i) {
    // ----------------------- 2016 -----------------------
    if(year==2016){
	    switch (mode){
	      case 0:
          Lumi_distr.push_back(Data2016Golden[i]);
          break;
	      case 1:
          Lumi_distr.push_back(Data2016Golden_up[i]);
	        break;
	      case -1:
          Lumi_distr.push_back(Data2016Golden_dn[i]);
	        break;
        default:
	        abort();
	        break;
	    } // end of switch
	    MC_distr.push_back(MC_2016_25ns_Moriond17[i]);
    }

    // ----------------------- 2017 -----------------------
    else if(year==2017){
      // to add 2017 PU distributions
	    switch (mode){
	      case 0:
          Lumi_distr.push_back(Data2016Golden[i]);
          break;
	      case 1:
          Lumi_distr.push_back(Data2016Golden_up[i]);
	        break;
	      case -1:
          Lumi_distr.push_back(Data2016Golden_dn[i]);
	        break;
        default:
	        abort();
	        break;
	    } // end of switch
	    MC_distr.push_back(MC_2016_25ns_Moriond17[i]);
    }
    else{
      std::cout << "Select a compatible year for PU reweighting." << std::endl;
    }
  } // end of loop over bins

  // no histograms for input: use vectors
  
  // now, make histograms out of them:

  // first, check they are the same size...

  if( MC_distr.size() != Lumi_distr.size() ){   
    std::cout << "MC_distr.size() = " << MC_distr.size() << std::endl;
    std::cout << "Lumi_distr.size() = " << Lumi_distr.size() << std::endl;
    std::cerr << "ERROR: standalone_LumiReWeighting: input vectors have different sizes. Quitting... \n";
  }


  weights_ = new TH1F(Form("luminumer_%d",mode),
    Form("luminumer_%d",mode),
    NBins,0., float(NBins));

  TH1F* den = new TH1F(Form("lumidenom_%d",mode),
    Form("lumidenom_%d",mode),
    NBins,0., float(NBins));

  for(int ibin = 1; ibin<NBins+1; ++ibin ) {
    weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
    den->SetBinContent(ibin,MC_distr[ibin-1]);
  }

  //std::cout << "Data Input " << std::endl;
  //for(int ibin = 1; ibin<NBins+1; ++ibin){
    //std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  //}
  //std::cout << "MC Input " << std::endl;
  //for(int ibin = 1; ibin<NBins+1; ++ibin){
    //std::cout << "   " << ibin-1 << " " << den->GetBinContent(ibin) << std::endl;
  //}

  // check integrals, make sure things are normalized

  float deltaH = weights_->Integral();
  if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
    weights_->Scale( 1.0/ weights_->Integral() );
  }
  float deltaMC = den->Integral();
  if(fabs(1.0 - deltaMC) > 0.02 ) {
    den->Scale(1.0/ den->Integral());
  }

  weights_->Divide( den );  // so now the average weight should be 1.0    

  double inte = 0;
  for(int ibin = 1; ibin < NBins+1; ++ibin){
    inte += weights_->GetBinContent(ibin) * MC_distr[ibin-1];
  }

  std::cout << "PU weight normalisation: " << inte << "\n";

}

standalone_LumiReWeighting::~standalone_LumiReWeighting()
{
}

double standalone_LumiReWeighting::weight( int npv ) {
  int bin = weights_->GetXaxis()->FindBin( npv );
  return weights_->GetBinContent( bin );
}

#endif
