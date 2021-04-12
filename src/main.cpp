#include <iostream>
#include <string>
#include "Kinematics.h"
#include "SabreEfficiency.h"
#include "SabreDetector.h"
#include "KinematicsExceptions.h"

int main(int argc, char** argv) {
	if(argc<2) {
		std::cerr<<"Incorrect number of arguments!"<<std::endl;
		return 1;
	}
	
	Mask::Kinematics calculator;
	try {
		if(!calculator.LoadConfig(argv[1])) {
			return 1;
		}
		calculator.Run();
	} catch(const std::exception& e) {
		std::cerr<<"Exception caught! Information: "<<e.what()<<std::endl;
		std::cerr<<"Terminating process."<<std::endl;
		return 1;
	}
	SabreEfficiency sabre;
	std::string mapfile = "./etc/DeadChannels.txt";
	sabre.SetReactionType(calculator.GetReactionType());
	sabre.SetDeadChannelMap(mapfile);
	sabre.CalculateEfficiency(calculator.GetOutputName());

	/*std::vector<SabreDetector> detectors;
	const double INNER_R = 0.0326;
    const double OUTER_R = 0.1351;
    const double TILT = 40.0;
    const double DIST_2_TARG = -0.1245;
    const double PHI_COVERAGE = 54.4; //delta phi for each det
    const double PHI0 = 234.0; //center phi values for each det in array
    const double PHI1 = 162.0; //# is equal to detID in channel map
    const double PHI2 = 306.0;
    const double PHI3 = 18.0;
    const double PHI4 = 90.0;
    const double DEG2RAD = M_PI/180.0;

    detectors.reserve(5);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI0*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI1*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI2*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI3*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI4*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);

 	double theta, phi, expected_flat_p;
 	for(int h=0; h<5; h++) {
 		for(int j=0; j<16; j++) {
 			for(int k=0; k<4; k ++) {
 				theta  = detectors[h].GetRingTiltCoords(j, k).GetTheta();
 				phi = detectors[h].GetRingTiltCoords(j, k).GetPhi();
 				expected_flat_p = detectors[h].GetRingFlatCoords(j, k).GetPhi();  
 				for(int i=0; i<5; i++) {
 					auto channels = detectors[i].GetTrajectoryRingWedge(theta, phi);
 					if(channels.first != -1) {
 						std::cout<<"Detected in detector"<<i<<" ring: "<<channels.first<<" wedge: "<<channels.second<<" Expected -- detector: "<<h<<" ring: "<<j;
 						if(k == 0 || k == 1) std::cout<<" wedge: 0"<<std::endl;
 						else std::cout<<" wedge: 7"<<std::endl;
 						break;
 					} else if(i == 4) {
 						std::cout<<" Not found! detector: "<<h<<" ring: "<<j<<" corner: "<<k<<" theta: "<<theta/DEG2RAD<<" phi: "<<phi/DEG2RAD<<" flat_p: "<<expected_flat_p/DEG2RAD<<std::endl;
 					}
 				}
 			}
 		}
 	}*/

	return 0;

}