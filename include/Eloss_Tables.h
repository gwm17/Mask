#ifndef ELOSS_TABLES_H
#define ELOSS_TABLES_H
	
namespace Mask {

	#define MAX_Z 93 //Maximum number of elements for which we have hydrogen coefficients

	/*Atomic Masses for elements H through U. Taken from ELAST data*/
	static double NATURAL_MASS[MAX_Z] =        
	                       {0,
	                        1.00797,   4.0026,   6.939,    9.0122,   10.818,
	                        12.01115,  14.0067,  15.9994,  18.99984,  20.183,
	                        22.9898,   24.312,   26.9815,  28.086,    30.9738,
	                        32.064,     35.453,  39.948,   39.102,    40.08, 
	                        44.956,     47.90,   50.942,   51.996,    54.938,
	                        55.847,     58.933,  58.71,    63.54,     65.37, 
	                        69.72,      72.59,   74.922,   78.96,     79.909,
	                        83.80,      85.47,   87.62,    88.909,    91.22, 
	                        92.906,     95.94,   98.,     101.07,    102.905,
	                        106.4,      107.87,  112.4,    114.82,    118.69,
	                        121.75,     127.60,  126.904,  131.3,     132.905,
	                        137.34,     138.91,  140.12,   140.907,   144.24,
	                        146.,       150.35,  151.96,   157.25,   158.924,
	                        162.50,     164.93,  167.26,   168.934,   173.04,
	                        174.97,     178.49,  180.948,  183.85,     186.2,
	                        190.2,      192.2,   195.09,   196.967,    200.59,
	                        204.37,     207.19,  208.98,   209.,       210.,
	                        222.,       223.,    226.,     227.,       232.038,
	                        231.,       238.03 };
	
	/*Stopping power coefficients for hydrogen into Z
	 *Taken from Ziegler, Hydrogen Stopping Powers and Ranges in All Elements, Volume 3 of
	 *The Stopping and Ranges of Ions in Matter, 1977, Pergamon Press Inc.
	 *Expanded table for method given in O.G. Spanc
	 *Note: some elements represent extrapolations when there was no data to fit
	 */
	static double HYDROGEN_COEFF[MAX_Z][12] = {
	  /*Blank*/{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
	  /*H*/{1.262,1.44,242.6,1.2E4,0.1159,0.0005099,5.436E4,-5.052,2.049,-0.3044,0.01966,-0.0004659},
	  /*He*/{1.229,1.397,484.5,5873,0.05225,0.00102,2.451E4,-2.158,0.8278,-0.1172,0.007259,-0.000166},
	  /*Li*/{1.411,1.6,725.6,3013,0.04578,0.00153,2.147E4,-0.5831,0.562,-0.1183,0.009298,-0.0002498},
	  /*Be*/{2.248,2.59,966,153.8,0.03475,0.002039,1.63E4,0.2779,0.1745,-0.05684,0.005155,-0.0001488},
	  /*B*/{2.474,2.815,1206,1060,0.02855,0.002549,1.345E4,-2.445,1.283,-0.2205,0.0156,-0.000393},
	  /*C*/{2.631,2.989,1445,957.2,0.02819,0.003059,1.322E4,-4.38,2.044,-0.3283,0.02221,-0.0005417},
	  /*N*/{2.954,3.35,1683,1900,0.02513,0.003569,1.179E4,-5.054,2.325,-0.3713,0.02506,-0.0006109},
	  /*O*/{2.652,3,1920,2000,0.0223,0.004079,1.046E4,-6.734,3.019,-0.4748,0.03171,-0.0007669},
	  /*F*/{2.085,2.352,2157,2634,0.01816,0.004589,8517,-5.571,2.449,-0.3781,0.02483,-0.0005919},
	  /*Ne*/{1.951,2.199,2393,2699,0.01568,0.005099,7353,-4.408,1.879,-0.2814,0.01796,-0.0004168},
	  /*Na*/{2.542,2.869,2628,1854,0.01472,0.005609,6905,-4.959,2.073,-0.3054,0.01921,-0.0004403},
	  /*Mg*/{3.792,4.293,2862,1009,0.01397,0.006118,6551,-5.51,2.266,-0.3295,0.02047,-0.0004637},
	  /*Al*/{4.154,4.739,2766,164.5,0.02023,0.006628,6309,-6.061,2.46,-0.3535,0.02173,-0.0004871},
	  /*Si*/{4.15,4.7,3329,550,0.01321,0.007138,6194,-6.294,2.538,-0.3628,0.0222,-0.0004956},
	  /*P*/{3.232,3.647,3561,1560,0.01267,0.007648,5942,-6.527,2.616,-0.3721,0.02267,-0.000504},
	  /*S*/{3.447,3.891,3792,1219,0.01211,0.008158,5678,-6.761,2.694,-0.3814,0.02314,-0.0005125},
	  /*Cl*/{5.047,5.714,4023,878.6,0.01178,0.008668,5524,-6.994,2.773,-0.3907,0.02361,-0.0005209},
	  /*Ar*/{5.731,6.5,4253,530,0.01123,0.009178,5268,-7.227,2.851,-0.4,0.02407,-0.0005294},
	  /*K*/{5.151,5.833,4482,545.7,0.01129,0.009687,5295,-7.44,2.923,-0.4094,0.02462,-0.0005411},
	  /*Ca*/{5.521,6.252,4710,553.3,0.01112,0.0102,5214,-7.653,2.995,-0.4187,0.02516,-0.0005529},
	  /*Sc*/{5.201,5.884,4938,560.9,0.009995,0.01071,4688,-8.012,3.123,-0.435,0.02605,-0.0005707},
	  /*Ti*/{4.862,5.496,5165,568.5,0.009474,0.01122,4443,-8.371,3.251,-0.4513,0.02694,-0.0005886},
	  /*V*/{4.48,5.055,5391,952.3,0.009117,0.01173,4276,-8.731,3.379,-0.4676,0.02783,-0.0006064},
	  /*Cr*/{3.983,4.489,5616,1336,0.008413,0.01224,3946,-9.09,3.507,-0.4838,0.02872,-0.0006243},
	  /*Mn*/{3.469,3.907,5725,1461,0.008829,0.01275,3785,-9.449,3.635,-0.5001,0.02961,-0.0006421},
	  /*Fe*/{3.519,3.963,6065,1243,0.007782,0.01326,3650,-9.809,3.763,-0.5164,0.0305,-0.00066},
	  /*Co*/{3.14,3.535,6288,1372,0.007361,0.01377,3453,-10.17,3.891,-0.5327,0.03139,-0.0006779},
	  /*Ni*/{3.553,4.004,6205,555.1,0.008763,0.01428,3297,-10.53,4.019,-0.549,0.03229,-0.0006957},
	  /*Cu*/{3.696,4.175,4673,387.8,0.02188,0.01479,3174,-11.18,4.252,-0.5791,0.03399,-0.0007314},
	  /*Zn*/{4.21,4.75,6953,295.2,0.006809,0.0153,3194,-11.57,4.394,-0.598,0.03506,-0.0007537},
	  /*Ga*/{5.041,5.697,7173,202.6,0.006725,0.01581,3154,-11.95,4.537,-0.6169,0.03613,-0.0007759},
	  /*Ge*/{5.554,6.3,6496,110,0.009689,0.01632,3097,-12.34,4.68,-0.6358,0.03721,-0.0007981},
	  /*As*/{5.323,6.012,7611,292.5,0.006447,0.01683,3024,-12.72,4.823,-0.6547,0.03828,-0.0008203},
	  /*Se*/{5.874,6.656,7395,117.5,0.007684,0.01734,3006,-13.11,4.965,-0.6735,0.03935,-0.0008425},
	  /*Br*/{5.611,6.335,8046,365.2,0.006244,0.01785,2928,-13.4,5.083,-0.6906,0.04042,-0.0008675},
	  /*Kr*/{6.411,7.25,8262,220,0.006087,0.01836,2855,-13.69,5.2,-0.7076,0.0415,-0.0008925},
	  /*Rb*/{5.694,6.429,8478,292.9,0.006087,0.01886,2855,-13.92,5.266,-0.714,0.04173,-0.0008943},
	  /*Sr*/{6.339,7.159,8693,330.3,0.006003,0.01937,2815,-14.14,5.331,-0.7205,0.04196,-0.0008962},
	  /*Y*/{6.407,7.234,8907,367.8,0.005889,0.01988,2762,-14.36,5.397,-0.7269,0.04219,-0.000898},
	  /*Zr*/{6.734,7.603,9120,405.2,0.005765,0.02039,2704,-14.59,5.463,-0.7333,0.04242,-0.0008998},
	  /*Nb*/{6.902,7.791,9333,442.7,0.005587,0.0209,2621,-16.22,6.094,-0.8225,0.04791,-0.001024},
	  /*Mo*/{6.425,7.248,9545,480.2,0.005367,0.02141,2517,-17.85,6.725,-0.9116,0.05339,-0.001148},
	  /*Tc*/{6.799,7.671,9756,517.6,0.005315,0.02192,2493,-17.96,6.752,-0.9135,0.05341,-0.001147},
	  /*Ru*/{6.108,6.887,9966,555.1,0.005151,0.02243,2416,-18.07,6.779,-0.9154,0.05342,-0.001145},
	  /*Rh*/{5.924,6.677,1.018E4,592.5,0.004919,0.02294,2307,-18.18,6.806,-0.9173,0.05343,-0.001143},
	  /*Pd*/{5.238,5.9,1.038E4,630,0.004758,0.02345,2231,-18.28,6.833,-0.9192,0.05345,-0.001142},
	  /*Ag*/{5.623,6.354,7160,337.6,0.01394,0.02396,2193,-18.39,6.86,-0.9211,0.05346,-0.00114},
	  /*Cd*/{5.814,6.554,1.08E4,355.5,0.004626,0.02447,2170,-18.62,6.915,-0.9243,0.0534,-0.001134},
	  /*In*/{6.23,7.024,1.101E4,370.9,0.00454,0.02498,2129,-18.85,6.969,-0.9275,0.05335,-0.001127},
	  /*Sn*/{6.41,7.227,1.121E4,386.4,0.004474,0.02549,2099,-19.07,7.024,-0.9308,0.05329,-0.001121},
	  /*Sb*/{7.5,8.48,8608,348,0.009074,0.026,2069,-19.57,7.225,-0.9603,0.05518,-0.001165},
	  /*Te*/{6.979,7.871,1.162E4,392.4,0.004402,0.02651,2065,-20.07,7.426,-0.9899,0.05707,-0.001209},
	  /*I*/{7.725,8.716,1.183E4,394.8,0.004376,0.02702,2052,-20.56,7.627,-1.019,0.05896,-0.001254},
	  /*Xe*/{8.231,9.289,1.203E4,397.3,0.004384,0.02753,2056,-21.06,7.828,-1.049,0.06085,-0.001298},
	  /*Cs*/{7.287,8.218,1.223E4,399.7,0.004447,0.02804,2086,-20.4,7.54,-1.004,0.05782,-0.001224},
	  /*Ba*/{7.899,8.911,1.243E4,402.1,0.004511,0.02855,2116,-19.74,7.252,-0.9588,0.05479,-0.001151},
	  /*La*/{8.041,9.071,1.263E4,404.5,0.00454,0.02906,2129,-19.08,6.964,-0.9136,0.05176,-0.001077},
	  /*Ce*/{7.489,8.444,1.283E4,406.9,0.00442,0.02957,2073,-18.43,6.677,-0.8684,0.04872,-0.001003},
	  /*Pr*/{7.291,8.219,1.303E4,409.3,0.004298,0.03008,2016,-17.77,6.389,-0.8233,0.04569,-0.0009292},
	  /*Nd*/{7.098,8,1.323E4,411.8,0.004182,0.03059,1962,-17.11,6.101,-0.7781,0.04266,-0.0008553},
	  /*Pm*/{6.91,7.786,1.343E4,414.2,0.00405,0.0311,1903,-16.45,5.813,-0.733,0.03963,-0.0007815},
	  /*Sm*/{6.728,7.58,1.362E4,416.6,0.003976,0.03161,1865,-15.79,5.526,-0.6878,0.0366,-0.0007077},
	  /*Eu*/{6.551,7.38,1.382E4,419,0.003877,0.03212,1819,-15.13,5.238,-0.6426,0.03357,-0.0006339},
	  /*Gd*/{6.739,7.592,1.402E4,421.4,0.003863,0.03263,1812,-14.47,4.95,-0.5975,0.03053,-0.0005601},
	  /*Tb*/{6.212,6.996,1.421E4,423.9,0.003725,0.03314,1747,-14.56,4.984,-0.6022,0.03082,-0.0005668},
	  /*Dy*/{5.517,6.21,1.44E4,426.3,0.003632,0.03365,1703,-14.65,5.018,-0.6069,0.03111,-0.0005734},
	  /*Ho*/{5.219,5.874,1.46E4,428.7,0.003498,0.03416,1640,-14.74,5.051,-0.6117,0.03141,-0.0005801},
	  /*Er*/{5.071,5.706,1.479E4,433,0.003405,0.03467,1597,-14.83,5.085,-0.6164,0.0317,-0.0005867},
	  /*Tm*/{4.926,5.542,1.498E4,433.5,0.003342,0.03518,1567,-14.91,5.119,-0.6211,0.03199,-0.0005933},
	  /*Yb*/{4.787, 5.386,1.517E4,435.9,0.003292,0.03569,1544,-15,5.153,-0.6258,0.03228,-0.0006},
	  /*Lu*/{4.893, 5.505,1.536E4,438.4,0.003243,0.0362,1521,-15.09,5.186,-0.6305,0.03257,-0.0006066},
	  /*Hf*/{5.028, 5.657,1.555E4,440.8,0.003195,0.03671,1499,-15.18,5.22,-0.6353,0.03286,-0.0006133},
	  /*Ta*/{4.738, 5.329,1.574E4,443.2,0.003186,0.03722,1494,-15.27,5.254,-0.64,0.03315,-0.0006199},
	  /*W*/{4.574, 5.144,1.593E4,442.4,0.003144,0.03773,1475,-15.67,5.392,-0.6577,0.03418,-0.0006426},
	  /*Re*/{5.2, 5.851,1.612E4,441.6,0.003122,0.03824,1464,-16.07,5.529,-0.6755,0.03521,-0.0006654},
	  /*Os*/{5.07, 5.704,1.63E4,440.9,0.003082,0.03875,1446,-16.47,5.667,-0.6932,0.03624,-0.0006881},
	  /*Ir*/{4.945, 5.563,1.649E4,440.1,0.002965,0.03926,1390,-16.88,5.804,-0.711,0.03727,-0.0007109},
	  /*Pt*/{4.476, 5.034,1.667E4,439.3,0.002871,0.03977,1347,-17.28,5.942,-0.7287,0.0383,-0.0007336},
	  /*Au*/{4.856, 5.46,1.832E4,438.5,0.002542,0.04028,1354,-17.02,5.846,-0.7149,0.0374,-0.0007114},
	  /*Hg*/{4.308, 4.843,1.704E4,487.8,0.002882,0.04079,1352,-17.84,6.183,-0.7659,0.04076,-0.0007925},
	  /*Tl*/{4.723, 5.311,1.722E4,537,0.002913,0.0413,1366,-18.66,6.52,-0.8169,0.04411,-0.0008737},
	  /*Pb*/{5.319, 5.982,1.74E4,586.3,0.002871,0.04181,1347,-19.48,6.857,-0.8678,0.04747,-0.0009548},
	  /*Bi*/{5.956, 6.7,1.78E4,677,0.00266,0.04232,1336,-19.55,6.871,-0.8686,0.04748,-0.0009544},
	  /*Po*/{6.158, 6.928,1.777E4,586.3,0.002812,0.04283,1319,-19.62,6.884,-0.8694,0.04748,-0.000954},
	  /*At*/{6.204, 6.979,1.795E4,586.3,0.002776,0.04334,1302,-19.69,6.898,-0.8702,0.04749,-0.0009536},
	  /*Rn*/{6.181, 6.954,1.812E4,586.3,0.002748,0.04385,1289,-19.76,6.912,-0.871,0.04749,-0.0009532},
	  /*Fr*/{6.949, 7.82,1.83E4,586.3,0.002737,0.04436,1284,-19.83,6.926,-0.8718,0.0475,-0.0009528},
	  /*Ra*/{7.506, 8.448,1.848E4,586.3,0.002727,0.04487,1279,-19.9,6.94,-0.8726,0.04751,-0.0009524},
	  /*Ac*/{7.649, 8.609,1.866E4,586.3,0.002697,0.04538,1265,-19.97,6.953,-0.8733,0.04751,-0.000952},
	  /*Th*/{7.71, 8.679,1.883E4,586.3,0.002641,0.04589,1239,-20.04,6.967,-0.8741,0.04752,-0.0009516},
	  /*Pa*/{7.407, 8.336,1.901E4,586.3,0.002603,0.0464,1221,-20.11,6.981,-0.8749,0.04752,-0.0009512},
	  /*U*/{7.29, 8.204,1.918E4,586.3,0.002573,0.04691,1207,-20.18,6.995,-0.8757,0.04753,-0.0009508}
	};
}
#endif
