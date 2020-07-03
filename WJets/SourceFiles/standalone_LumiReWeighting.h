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
// Here we put the arrays of data points that represent the pileup profiles in data and in the corresponding MC
// The script will then normalize each the pileup profiles from data and MC, and then divide data/MC to get the weights dist.
// Systematic variations come from uncertainty in the pp min bias xsec used to calculate avg. number of pileup vertices per lumi section
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
double Data2016Golden_up[100] = {
  232682.41284357742, 659475.2271608508, 2183864.732392254, 2745785.9351763777, 4070508.568027147, 5392709.986098041, 6381715.140083467, 9038775.190907836, 23775131.93747618, 53984515.58845366, 
  115978344.95937279, 245979886.66618606, 443258734.2732075, 679924868.931169, 937608397.69176, 1218455780.0632157, 1501218416.8142667, 1729889027.479183, 1882024464.658047, 1968112810.7867124, 
  2012219732.5181403, 2040066909.5023842, 2054543180.8271363, 2036793049.4297233, 1978244564.1186273, 1889392659.616465, 1781221832.266064, 1657923662.5715144, 1523121501.0582304, 1382226716.120118, 
  1239652247.3950613, 1097651524.711982, 957682414.192569, 821902874.4640143, 693281727.1996627, 574788841.5294858, 468602491.9800088, 375738277.1572735, 296147312.2512948, 229053337.3642116, 
  173363343.06884006, 127950797.6314024, 91747750.52470127, 63700453.91348277, 42700047.358398184, 27570598.13794411, 17117074.677403446, 10204881.017004691, 5837000.179664271, 3201533.5920285764, 
  1683863.6674764492, 849882.6186800991, 412537.27082202036, 193593.24290710725, 88881.05371006558, 40970.17923951562, 19934.563272018906, 11012.539026035098, 7298.601154429589, 5722.910664296525, 
  4985.457537627968, 4560.493605597344, 4244.753362476094, 3963.5310129579475, 3691.484843073756, 3421.626370304515, 3153.55690616613, 2889.087616194446, 2630.6133851676886, 2380.5008889081087, 
  2140.85313889621, 1913.4192030060615, 1699.5648488509023, 1500.2725960812518, 1316.1585167737028, 1147.4998487914347, 994.2699030090865, 856.1776652473462, 732.7099258795139, 623.1740669248635, 
  526.7399124607468, 442.47933245231235, 369.4025817760287, 306.4906452815245, 252.7231344004297, 207.10153068012332, 168.66778871978838, 136.5184902482223, 109.8148802392202, 87.78921528673617, 
  69.74791634465632, 55.072046348367344, 43.215633228745396, 33.70233601304743, 26.120911818050153, 20.119890066492008, 15.401802169541996, 11.717254477249574, 8.859072966097328, 6.656692530037995
};

// Systematic variation down -- 69.2 mb * (1-0.046)
double Data2016Golden_down[100] = {
  247406.86190125695, 1069192.5329278684, 2428216.350572554, 3566711.2688976275, 4987870.597469167, 6585313.369619577, 8093967.8735435065, 19953541.876027748, 51890982.76976208, 119725462.67371367, 
  272658060.8250202, 512969483.011789, 802044185.1492237, 1118356640.9746308, 1462439983.1099296, 1779534467.9739075, 2004432905.4475088, 2134348205.2395785, 2197908047.8114195, 2233560909.039848, 
  2252302228.753103, 2229565382.0692945, 2154133531.3249874, 2041653621.406011, 1905607141.8185542, 1751102360.4194503, 1584625258.7182012, 1413637828.308284, 1242566424.6290426, 1073674033.3579733, 
  909807232.7607071, 755134318.9525062, 613914677.0516247, 489145748.0897949, 382000174.4776026, 292073214.82415336, 218027697.99091268, 158235274.30351537, 111116205.82766648, 75149256.12556043, 
  48754819.66959769, 30247030.916957725, 17901664.75746416, 10090410.045706483, 5410591.292645974, 2758516.970862317, 1337572.920506894, 617864.2626171438, 273170.8994969539, 116956.92524208047, 
  49866.176360217825, 22468.145526396303, 11744.542953377517, 7639.007197198009, 6018.734040666122, 5281.033408992601, 4837.60457449107, 4484.04093876899, 4154.076913217944, 3829.091369189568, 
  3506.2677569396365, 3187.727822851842, 2876.8648763339424, 2577.0938946067313, 2291.414909031725, 2022.2611079082312, 1771.4569554606624, 1540.2266982422593, 1329.2308798916517, 1138.62052717332, 
  968.1026370492394, 817.0120688464161, 684.3857268190642, 569.0355597440423, 469.61755049579415, 384.6945338856709, 312.7913410764457, 252.4413919945654, 202.22441541526564, 160.79544939742547, 
  126.90565069508874, 99.41571710299677, 77.30290495136639, 59.662714258298294, 45.706329312250936, 34.75485777663427, 26.23132248899052, 19.651242061265236, 14.612502714336427, 10.785085902355434, 
  7.901083151729225, 5.745307630266487, 4.146705419734755, 2.9706804069739423, 2.112375598482057, 1.4908997092762233, 1.0444494273398952, 0.7262526435757266, 0.501243776481326, 0.3433767281709743
};

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
double Data2017Golden_up[100] = {
  250717.5593711223, 936720.9855645779, 1647697.344123416, 3798693.261071138, 3493737.093451203, 5208332.526360814, 6232770.150256079, 6253652.232822476, 7454669.239129748, 14162020.671089714, 
  31169279.99526131, 57552519.19691123, 99875033.26659426, 146015729.57951295, 204475579.81175125, 283995760.84713316, 393750627.7141839, 537633143.7857208, 698695818.2410258, 855570494.134949, 
  1000036948.2495108, 1126819196.1728737, 1229806038.6850753, 1309509186.7727358, 1371947836.2002254, 1425338988.950387, 1475957023.809288, 1521087938.625004, 1552557445.0446072, 1564082356.0098047, 
  1553037334.4808366, 1520471012.0937393, 1469982559.5188708, 1404850427.6218965, 1326957791.3409374, 1238951634.0234797, 1145722509.6534429, 1052949511.3008932, 964867183.8964175, 883875518.4320136, 
  811869974.0411644, 751415090.4854723, 705624864.9742507, 677174810.185184, 667256448.5936446, 674865788.7310822, 696340423.5152595, 725220219.104561, 752661324.1937772, 768765543.3792934, 
  764727265.3557508, 735070206.0839245, 678934530.1367173, 600903099.512204, 509164718.2737001, 413253406.6028359, 321892540.42033976, 241344913.19131994, 174828426.05438465, 122866469.71322331, 
  84134120.55047622, 56372434.85957338, 37106579.88312415, 24083431.78138701, 15464239.791301757, 9854695.12554707, 6251354.291924594, 3959304.7556308955, 2511093.7834888995, 1599323.9082774064, 
  1025489.3084203491, 663278.175740401, 433254.0340748889, 285890.06197837647, 190467.33381706456, 127956.58044861288, 86534.2258545219, 58799.86376647015, 40069.679938606474, 27337.700223067626, 
  18645.15343312981, 12696.389210617785, 8622.816268601557, 5835.707498968175, 3932.731276732553, 2637.412072428393, 1759.1619926705957, 1166.4488539585616, 768.5411004142684, 502.9652487560532, 
  326.8306995199296, 210.8052937144313, 134.9236913561036, 85.67038642427515, 53.95223863061777, 33.692755503855004, 20.86097508314, 12.803712463964358, 7.788989556660058, 4.695894179875592
};

// Systematic variation down -- 69.2 mb * (1-0.046)
double Data2017Golden_down[100] = {
  270631.23324132094, 1199612.637835579, 2624462.286788161, 3652948.322194754, 4674229.912228364, 6554084.806432281, 6729741.9650886785, 7654125.105669664, 13360679.257476006, 32430321.503089454, 
  63798464.356859334, 114824692.03345194, 171881960.5554621, 247898324.60218096, 355057093.13371897, 505960358.05305904, 693414704.317184, 885468631.2763422, 1063487310.3690034, 1219911178.5643892, 
  1345864785.447774, 1441187565.7516923, 1514539008.0023394, 1577759402.7363992, 1637290010.017443, 1685812118.6786244, 1712058464.166281, 1709623798.9592476, 1677929414.5797048, 1621315213.3278198, 
  1544811794.7964318, 1451373702.8010893, 1344879764.8318093, 1232534784.9523563, 1122250337.215936, 1019333733.8649449, 926787056.3404257, 847731657.7965287, 786414951.6757892, 747046370.9856322, 
  731933897.6899104, 740104411.3597095, 766443519.109516, 801439408.6593196, 832005942.8056625, 844018309.2695035, 826185008.1059879, 773332933.9391261, 688469733.9914783, 581788746.4066515, 
  466740496.8670874, 356262099.200293, 259717890.0678228, 181715888.72106066, 122694173.5626463, 80399531.05538642, 51412213.62653704, 32245905.560298517, 19929123.165051553, 12188464.690929625, 
  7406310.419049323, 4489102.486743885, 2724696.4678632407, 1662261.3439648533, 1022634.4894851497, 635961.9815010906, 400291.51340258843, 255002.1736102726, 164217.89328921065, 106692.68134850933, 
  69764.17955577976, 45796.750068990135, 30112.25751823119, 19792.25099686113, 12982.918843614505, 8487.656911263768, 5524.092058160846, 3575.952746116776, 2300.5955596709655, 1469.9772210547324, 
  932.2687300031654, 586.5408000073166, 365.9087968416385, 226.24463375165698, 138.59536765941553, 84.08873634515857, 50.514654012547055, 30.038319783055734, 17.677288284957516, 10.293274408444045, 
  5.929483976201004, 3.378648749497203, 1.904047080071352, 1.0611421180121892, 0.584775649419019, 0.3186312656689829, 0.1716478855424557, 0.09141351615178911, 0.04812579665491834, 0.025044765143582004
};

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// For 2018 data, All Eras --------------------------
// entries go from 0 vertices to 99 vertices

// Central value -- 69.2 mb
double Data2018Golden[100] = {
  271646.1851691027, 956235.9896633611, 3039372.3064388265, 6828821.712880966, 12155546.162411349, 18945105.709091246, 28001037.602586146, 40768339.75487006, 57722392.70738312, 80314046.68547472, 
  112501851.61579308, 158170667.54948705, 220442568.7900867, 301591592.13116646, 403483910.42947644, 527959448.90674955, 675904258.9329485, 845775594.1396854, 1032276549.9556235, 1225851909.628376, 
  1413734557.5090747, 1583029703.596016, 1724838893.8850346, 1837112512.1903586, 1924173706.1644573, 1993601915.1701465, 2052969165.5496466, 2107816428.0124907, 2160932250.817238, 2212635728.663044, 
  2261421417.759793, 2304731104.6378593, 2339534441.6951456, 2362572706.5900187, 2370325223.704171, 2359037548.906971, 2324868902.477417, 2264395535.687464, 2175293581.5565586, 2057052401.2082696, 
  1911459896.376943, 1742695571.5987933, 1557004932.4078088, 1362019301.0410178, 1165872659.3414512, 976325668.0544866, 799962154.882498, 641642305.3318702, 504264138.47279644, 388793557.47466624, 
  294559452.9744404, 219701513.14549437, 161651994.96749765, 117573936.67569274, 84694339.07977186, 60522242.455085516, 42954016.99511572, 30296301.310013723, 21237164.28150894, 14788592.143966574, 
  10221544.761802876, 7004973.7689238675, 4754612.315404755, 3193011.559912291, 2119854.392708486, 1390551.1558580967, 900989.5714379459, 576629.8041423364, 364595.78066639154, 227840.92483313987, 
  140789.87478025994, 86071.23072749135, 52083.496230155106, 31207.541019569628, 18519.525277135002, 10884.889281279484, 6335.367588475792, 3650.279903151363, 2081.0203569287783, 1173.1812858435978, 
  653.5927179752027, 359.58862382834997, 195.23940869186325, 104.5469022357118, 55.179331061459315, 28.68973799209806, 14.687578255012443, 7.400495037966902, 3.668557720342649, 1.7885972855044787, 
  0.8574138963987736, 0.40404057996929593, 0.18712349620570948, 0.08515888050590181, 0.0380785229179208, 0.016728652496056796, 0.007221276943880962, 0.00306411722748198, 0.0012793215725774927, 0.0005268888877819574
};

// Systematic variation up -- 69.2 mb * (1+0.046)
double Data2018Golden_up[100] = {
  252527.3942835556, 835039.9531189795, 2595656.676559611, 5868311.017598032, 10486216.951094942, 16409286.146796683, 23960798.669375796, 34500698.34686839, 48616350.38773735, 66866154.35970679, 
  91902419.83618587, 127275095.44546, 175846895.5277399, 239930620.910274, 321166051.5435692, 421137659.80002403, 541249784.6676397, 681793645.7632687, 840791619.3760328, 1013080941.6213002, 
  1190120997.690366, 1361021922.250256, 1515094771.4674127, 1645029946.4551826, 1749052167.5542479, 1830577016.1062436, 1895839248.5990827, 1951325818.9546833, 2002069530.4118357, 2050931420.2381854, 
  2098689786.5559638, 2144498632.449478, 2186456808.1313195, 2222153306.340235, 2248942438.0840616, 2263997689.3635755, 2264280581.520892, 2246592735.0737066, 2207784482.5534554, 2145224204.482576, 
  2057335311.3424034, 1944124385.0298643, 1807483216.2532835, 1651198280.2763958, 1480653947.497025, 1302286533.6888385, 1122902205.7049468, 949009285.6415987, 786197172.9720702, 638703016.5651118, 
  509205181.4536392, 398815323.1668486, 307267724.518257, 233245947.07848716, 174751121.0516691, 129455446.27260706, 94988843.198578, 69143451.34767698, 49991416.9282622, 35930194.82901639, 
  25679701.57904631, 18248734.476825126, 12887289.361964228, 9037076.825792773, 6286656.675555232, 4334290.0815860415, 2958997.9761201944, 1998926.739623989, 1335569.4988241557, 882362.3124604418, 
  576406.7500806799, 372381.63044861215, 237993.68885974967, 150537.0896708219, 94280.41835125738, 58490.95119961222, 35958.45810274715, 21911.0107825364, 13234.630788425007, 7923.528865744283, 
  4700.939106671165, 2762.8179714883818, 1607.7453855274698, 925.8549128121418, 527.3204215828362, 296.8604925857506, 165.09036658103847, 90.644339006324, 49.11170744427789, 26.245256789154592, 
  13.82786980416477, 7.180197834919491, 3.673261620367264, 1.8508761160033815, 0.9183431519830497, 0.4485819478156339, 0.21567939896755364, 0.10205646275367414, 0.04752131250957054, 0.02177341074692686
};

// Systematic variation down -- 69.2 mb * (1-0.046)
double Data2018Golden_down[100] = {
  293878.8351269457, 1105698.5885746558, 3583797.961871093, 7990315.87212832, 14171961.370980725, 22073185.628542222, 33103851.8488191, 48655742.10350035, 69327201.86446178, 98168022.50268891, 
  140241478.8150954, 199742184.3819281, 279931331.65221006, 383148939.11479294, 511678208.99848413, 667115880.4539194, 848563909.1215931, 1050750648.5598707, 1263051979.378605, 1470464952.4567742, 
  1657329541.054008, 1812716760.0792823, 1934249834.797748, 2027439155.7434983, 2101542902.4940865, 2165397944.7426424, 2225025787.004445, 2282984759.853094, 2338919904.098592, 2390476246.64244, 
  2434290844.9602175, 2466570199.4221263, 2483230782.8911767, 2479859829.754949, 2451824681.227357, 2394732113.4527464, 2305282647.103318, 2182218169.107687, 2027074739.4998782, 1844422601.416359, 
  1641525992.645846, 1427494872.252796, 1212129810.1965835, 1004758888.9546161, 813185596.8289143, 642994013.8903868, 497277266.48267007, 376738076.53082806, 280137187.65581304, 204907913.98673558, 
  147785404.8929442, 105341007.88497768, 74362638.47707152, 52072725.26719611, 36207711.76096404, 25007261.456423543, 17149945.777324535, 11668711.990447603, 7867543.22464715, 5249971.648454767, 
  3463037.8334131767, 2255886.066187808, 1450262.371105237, 919816.4513431343, 575541.2689067504, 355375.1378604641, 216638.18323037893, 130458.29308590866, 77652.687603447, 45710.63595446466, 
  26620.537541653957, 15340.142853803205, 8746.405010696297, 4932.804290683357, 2750.5098954819173, 1515.348594529522, 824.2782269319811, 442.34203366729525, 234.00503662023047, 121.94029002863688, 
  62.54885424032836, 31.562028329603017, 15.658044421952319, 7.63345252987321, 3.655338649938383, 1.7186822132261217, 0.7932069259358306, 0.3592386734019648, 0.1596202657019684, 0.06957051574174697, 
  0.02974036335084787, 0.012469500079482923, 0.005128909148577475, 0.0020709793014500027, 0.0008224275782806997, 0.0003226814183310125, 0.00012646932755744822, 5.077548777554449e-05, 2.1958926201832815e-05, 1.10149865319059e-05
};

// For 2018 MC: RunIIAutumn18MiniAOD,
// pileup profile used is 102X_upgrade2018_realistic_v15
double MC_2018_RunIIAutumn18MiniAOD[100] = {
  0.0, 962.194504528577, 664.4551472304223, 4715.803093745268, 11866.931192945201, 25245.211863038396, 50376.15481233017, 97321.58462427967, 182444.38358802354, 337903.61110599746, 
  584666.5572911192, 954642.9057275312, 1485298.7059914763, 2203761.0727361096, 3078818.285537084, 4148577.1614175425, 5368877.148517095, 6724587.884063044, 8121360.442732052, 9624643.080032095, 
  11151329.07147432, 12671925.062503127, 14168893.469025912, 15620454.994026262, 17048555.149213973, 18365904.75081106, 19583202.704973586, 20774342.63238058, 21843163.07912159, 22872696.99662313, 
  23667033.042789347, 24401312.855962563, 24915436.779763795, 25262410.470749628, 25413473.804182157, 25424294.173340887, 25176597.90244323, 24775279.76290737, 24197155.756107442, 23428372.151275687, 
  22475939.38329499, 21427011.872771688, 20282327.02648749, 19018580.11837783, 17678052.228460617, 16244327.522242634, 14861766.856876943, 13441035.866472743, 12029610.531040706, 10728431.414778113, 
  9402506.85457499, 8154690.720953081, 7033290.537785603, 5955717.41818012, 5008153.567918281, 4163065.3276565014, 3421217.3318768465, 2794824.3267669035, 2249735.927044986, 1782980.2034363654, 
  1404275.164318329, 1101684.1935984057, 859185.4094488326, 661233.4729341257, 505815.93021842453, 382390.4235341235, 290776.6016572521, 215577.95479123617, 160045.74397605544, 117474.0795949867, 
  85551.9033347071, 61367.276247586175, 44784.28179917841, 33270.95626033091, 21743.814093632598, 16338.443905642456, 10688.862872919646, 7000.418294376378, 5466.656805222003, 2893.9534261583876, 
  2535.4536857990547, 1613.113639631185, 1110.7597065762047, 656.8511220209683, 464.2851786784007, 283.02722450914536, 322.4252503181061, 147.44025509565097, 186.7255884914979, 173.8105576851532, 
  89.97460549235446, 90.038016706002, 51.52008555906235, 12.757888367408459, 25.440131096914605, 6.4924126405576965, 25.684397141628864, 38.8378428295562, -0.07564563790231116, 51.59573119696466
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

standalone_LumiReWeighting::standalone_LumiReWeighting(int year, int mode){
  
  std::vector<float> Lumi_distr;
  std::vector<float> MC_distr;
  Lumi_distr.clear();
  MC_distr.clear();

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
      std::cout << "Using +1 sigma 4.6% value " << std::endl;
      break;
    case -1:
      std::cout << "Using -1 sigma 4.6% value " << std::endl;
      break;
    default:
      std::cout << "Using central value " << std::endl;
      break;
  } // end of switch

  Int_t NBins = 100;
  for (int i(0); i < NBins; ++i){
    // ----------------------- 2016 -----------------------
    if (year == 2016){
	    switch (mode){
	      case 0:
          Lumi_distr.push_back(Data2016Golden[i]);
          break;
	      case 1:
          Lumi_distr.push_back(Data2016Golden_up[i]);
	        break;
	      case -1:
          Lumi_distr.push_back(Data2016Golden_down[i]);
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
          Lumi_distr.push_back(Data2017Golden_up[i]);
	        break;
	      case -1:
          Lumi_distr.push_back(Data2017Golden_down[i]);
	        break;
        default:
	        abort();
	        break;
	    } // end of switch
	    MC_distr.push_back(MC_2017_RunIIFall17MiniAODv2[i]);
    }
    // ----------------------- 2018 -----------------------
    else if (year == 2018){
      switch (mode){
	      case 0:
          Lumi_distr.push_back(Data2018Golden[i]);
          break;
	      case 1:
          Lumi_distr.push_back(Data2018Golden_up[i]);
	        break;
	      case -1:
          Lumi_distr.push_back(Data2018Golden_down[i]);
	        break;
        default:
	        abort();
	        break;
	    } // end of switch
	    MC_distr.push_back(MC_2018_RunIIAutumn18MiniAOD[i]);
    }
    else std::cout << "Select a compatible year for PU reweighting." << std::endl;
  } // end of loop over bins

  // Check that vectors are of same length
  if ( MC_distr.size() != Lumi_distr.size() ){   
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
  //     std::cout << "   # vtx = " << ibin-1 << ": " << weights_->GetBinContent(ibin) << std::endl;
  //     // std::cout << weights_->GetBinContent(ibin) << std::endl;
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
