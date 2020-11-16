#include <iostream>
#include "Kinematics.h"
#include "SabreEfficiency.h"
#include "KinematicsExceptions.h"

#include <TGraph2D.h>
#include <TGraph.h>
#include <TFile.h>
#include <TAxis.h>
#include <TCanvas.h>

int main(int argc, char** argv) {
	if(argc<2) {
		std::cerr<<"Incorrect number of arguments!"<<std::endl;
		return 1;
	}
	
	/*Kinematics calculator;
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
	sabre.SetReactionType(calculator.GetReactionType());
	sabre.CalculateEfficiency(calculator.GetOutputName());*/

	std::vector<SabreDetGeometry> detectors;
	const double INNER_R = 0.0326;
    const double OUTER_R = 0.1351;
    const double TILT = 40.0;
    const double DIST_2_TARG = 0.1245;
    const double PHI_COVERAGE = 54.4; //delta phi for each det
    const double PHI0 = 36.0; //center phi values for each det in array
    const double PHI1 = 108.0; //# is equal to detID in channel map
    const double PHI2 = 324.0;
    const double PHI3 = 252.0;
    const double PHI4 = 180.0;
    const double DEG2RAD = M_PI/180.0;

    detectors.reserve(5);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI0*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI1*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI2*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI3*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI4*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);

 	double theta = 145*M_PI/180.0; double phi = 92.1428*M_PI/180.0;
 	std::cout<<"theta: "<<theta<<" phi: "<<phi<<std::endl;
 	for(int i=0; i<5; i++) {
 		if(detectors[i].IsInside(theta, phi)) {
 			std::cout<<"Caught it in detector: "<<i<<std::endl;
 		}
 		break;
 	}


	return 0;

}