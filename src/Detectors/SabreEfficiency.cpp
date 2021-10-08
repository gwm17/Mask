#include "SabreEfficiency.h"
#include "MaskFile.h"
#include <fstream>
#include <iostream>
#include <iomanip>


SabreEfficiency::SabreEfficiency() : 
	DetectorEfficiency(), deadlayer(DEADLAYER_THIN), sabre_eloss(SABRE_THICKNESS)
{
	detectors.reserve(5);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI0*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI1*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI2*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI3*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);
 	detectors.emplace_back(INNER_R,OUTER_R,PHI_COVERAGE*DEG2RAD,PHI4*DEG2RAD,TILT*DEG2RAD,DIST_2_TARG);

 	
	std::vector<int> dead_z = {14};
	std::vector<int> dead_a = {28};
	std::vector<int> dead_stoich = {1};
	deadlayer.SetElements(dead_z, dead_a, dead_stoich);
	sabre_eloss.SetElements(dead_z, dead_a, dead_stoich);
}

SabreEfficiency::~SabreEfficiency() {}

void SabreEfficiency::CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) {
	std::cout<<"----------SABRE Efficiency Calculation----------"<<std::endl;
	std::cout<<"Loading in output from kinematics simulation: "<<inputname<<std::endl;
	std::cout<<"Running efficiency calculation..."<<std::endl;

	if(!dmap.IsValid()) {
		std::cerr<<"Unable to run SABRE Efficiency without a dead channel map."<<std::endl;
		std::cerr<<"If you have no dead channels, simply make a file that's empty"<<std::endl;
		std::cerr<<"Exiting."<<std::endl;
		std::cout<<"---------------------------------------------"<<std::endl;
	}


	Mask::MaskFile input(inputname, Mask::MaskFile::FileType::read);
	Mask::MaskFile output(outputname, Mask::MaskFile::FileType::read);
	std::ofstream stats(statsname);
	stats<<std::setprecision(5);

	Mask::MaskFileHeader header = input.ReadHeader();
	output.WriteHeader(header.rxn_type, header.nsamples);
	stats<<"Efficiency statistics for data from "<<inputname<<" using the SPS-SABRE geometry"<<std::endl;
	stats<<"Given in order of target=0, projectile=1, ejectile=2, residual=3, .... etc."<<std::endl;

	Mask::MaskFileData data;
	Mask::Nucleus nucleus;

	std::vector<int> counts;
	switch(header.rxn_type) {
		case Mask::RxnType::PureDecay:
			counts.resize(3, 0);
			break;
		case Mask::RxnType::OneStepRxn:
			counts.resize(4, 0);
			break;
		case Mask::RxnType::TwoStepRxn:
			counts.resize(6, 0);
			break;
		case Mask::RxnType::ThreeStepRxn:
			counts.resize(8, 0);
			break;
		default:
		{
			std::cerr<<"Bad reaction type at AnasenEfficiency::CalculateEfficiency (given value: "<<GetStringFromRxnType(header.rxn_type)<<"). Quiting..."<<std::endl;
			input.Close();
			output.Close();
			stats.close();
			return;
		}
	}

	int percent5 = header.nsamples*0.05;
	int count = 0;
	int npercent = 0;

	while(true) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		data = input.ReadData();
		if(data.eof)
			break;

		for(unsigned int i=0; i<data.Z.size(); i++) {
			nucleus.SetIsotope(data.Z[i], data.A[i]);
			nucleus.SetVectorSpherical(data.theta[i], data.phi[i], data.p[i], data.E[i]);
			auto result = IsSabre(nucleus);
			if(result.first) {
				data.detect_flag[i] = true;
				counts[i]++;
			} else if(data.detect_flag[i] == true) {
				data.detect_flag[i] = false;
			}
		}

		output.WriteData(data);
	}

	input.Close();
	output.Close();

	stats<<std::setw(10)<<"Index"<<"|"<<std::setw(10)<<"Efficiency"<<std::endl;
	stats<<"---------------------"<<std::endl;
	for(unsigned int i=0; i<counts.size(); i++) {
		stats<<std::setw(10)<<i<<"|"<<std::setw(10)<<((double)counts[i]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	stats.close();

	std::cout<<std::endl;
	std::cout<<"Complete."<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
}

void SabreEfficiency::DrawDetectorSystem(const std::string& filename) {
	std::ofstream output(filename);

	std::vector<double> ringxs, ringys, ringzs;
    std::vector<double> wedgexs, wedgeys, wedgezs;
	Mask::Vec3 coords;
 	for(int i=0; i<5; i++) {
 		for(int j=0; j<detectors[i].GetNumberOfRings(); j++) {
 			for(int k=0; k<4; k++) {
 				coords = detectors[i].GetRingTiltCoords(j, k);
 				ringxs.push_back(coords.GetX());
 				ringys.push_back(coords.GetY());
 				ringzs.push_back(coords.GetZ());
 			}
 		}
 	}

 	for(int i=0; i<5; i++) {
 		for(int j=0; j<detectors[i].GetNumberOfWedges(); j++) {
 			for(int k=0; k<4; k++) {
 				coords = detectors[i].GetWedgeTiltCoords(j, k);
 				wedgexs.push_back(coords.GetX());
 				wedgeys.push_back(coords.GetY());
 				wedgezs.push_back(coords.GetZ());
 			}
 		}
 	}

 	output<<"SABRE Geometry File -- Coordinates for Detectors"<<std::endl;
	output<<"Edges: x y z"<<std::endl;
	for(unsigned int i=0; i<ringxs.size(); i++) {
		output<<ringxs[i]<<" "<<ringys[i]<<" "<<ringzs[i]<<std::endl;
	}
	for(unsigned int i=0; i<wedgexs.size(); i++) {
		output<<wedgexs[i]<<" "<<wedgeys[i]<<" "<<wedgezs[i]<<std::endl;
	}

 	output.close();
}

double SabreEfficiency::RunConsistencyCheck() {
	double theta, phi;
	double npoints = 5.0*16.0*4.0;
	int count=0;
 	for(int h=0; h<5; h++) {
 		for(int j=0; j<16; j++) {
 			for(int k=0; k<4; k ++) {
 				theta  = detectors[h].GetRingTiltCoords(j, k).GetTheta();
 				phi = detectors[h].GetRingTiltCoords(j, k).GetPhi();
 				for(int i=0; i<5; i++) {
 					auto channels = detectors[i].GetTrajectoryRingWedge(theta, phi);
 					if(channels.first != -1) {
 						count++;
 					}
 				}
 			}
 		}
 	}

 	return ((double)count)/npoints;
}

/*Returns if detected, as well as total energy deposited in SABRE*/
std::pair<bool,double> SabreEfficiency::IsSabre(Mask::Nucleus& nucleus) {
	if(nucleus.GetKE() <= ENERGY_THRESHOLD) {
		return std::make_pair(false, 0.0);
	}

	Mask::Vec3 coords;
	double thetaIncident, eloss, e_deposited;
	for(int i=0; i<5; i++) {
		auto chan = detectors[i].GetTrajectoryRingWedge(nucleus.GetTheta(), nucleus.GetPhi());
		if(chan.first != -1 && chan.second != -1) {
			if(dmap.IsDead(i, chan.first, 0) || dmap.IsDead(i, chan.second, 1)) break; //dead channel check
			coords = detectors[i].GetTrajectoryCoordinates(nucleus.GetTheta(), nucleus.GetPhi());
			thetaIncident = std::acos(coords.Dot(detectors[i].GetNormTilted())/(coords.GetR()));
			eloss = deadlayer.GetEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), M_PI - thetaIncident);
			if((nucleus.GetKE() - eloss) <= ENERGY_THRESHOLD) break; //deadlayer check
			e_deposited = sabre_eloss.GetEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE() - eloss, M_PI - thetaIncident);
			return std::make_pair(true, e_deposited);
		}
	}

	return std::make_pair(false,0.0);
}
