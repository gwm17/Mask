#include <iostream>
#include <string>
#include "Kinematics.h"
#include "SabreEfficiency.h"
#include "AnasenEfficiency.h"
#include "KinematicsExceptions.h"

#include <TGraph2D.h>
#include <TCanvas.h>

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
	
	

	/*
	SabreEfficiency sabre;
	std::string mapfile = "./etc/DeadChannels.txt";
	sabre.SetReactionType(calculator.GetReactionType());
	sabre.SetDeadChannelMap(mapfile);
	sabre.CalculateEfficiency(calculator.GetOutputName());
	//std::cout<<"Running consistency check(1=success): "<<sabre.RunConsistencyCheck()<<std::endl;
	//sabre.DrawDetectorSystem("/data1/gwm17/10B3He/Feb2021/simulation/SABREGeo.root");
	*/


 	AnasenEfficiency anasen;
 	anasen.SetReactionType(calculator.GetReactionType());
 	anasen.CalculateEfficiency(calculator.GetOutputName());
 	//std::cout<<"Running consistency check(1=success): "<<anasen.RunConsistencyCheck()<<std::endl;
 	//anasen.DrawDetectorSystem("/data1/gwm17/TRIUMF_7Bed/simulation/ANASENGeo_centered_target_targetGap_BackQQQ_test.root");


	return 0;

}