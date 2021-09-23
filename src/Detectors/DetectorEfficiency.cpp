#include "SabreEfficiency.h"
#include "AnasenEfficiency.h"
#include "KinematicsExceptions.h"
#include <iostream>
#include <string>

int main(int argc, char** argv) {

	if(argc != 4) {
		std::cerr<<"Incorrect number of commandline arguments! Returning."<<std::endl;
		return 1;
	}

	std::string inputname = argv[1];
	std::string outputname = argv[2];
	std::string statsname = argv[3];

	/*
	SabreEfficiency sabre;
	std::string mapfile = "./etc/DeadChannels.txt";
	sabre.SetDeadChannelMap(mapfile);
	sabre.CalculateEfficiency(inputname, outputname, statsname);
	//std::cout<<"Running consistency check(1=success): "<<sabre.RunConsistencyCheck()<<std::endl;
	//sabre.DrawDetectorSystem("/data1/gwm17/10B3He/Feb2021/simulation/SABREGeo.txt");
	*/

	try {
 	AnasenEfficiency anasen;
 	anasen.CalculateEfficiency(inputname, outputname, statsname);
 	//std::cout<<"Running consistency check(1=success): "<<anasen.RunConsistencyCheck()<<std::endl;
 	//anasen.DrawDetectorSystem("/data1/gwm17/TRIUMF_7Bed/simulation/ANASENGeo_centered_target_targetGap_BackQQQ_test.root");
	} catch(const std::exception& e) {
		std::cerr<<"Error: "<<e.what()<<std::endl;
		std::cerr<<"Terminating."<<std::endl;
		return 1;
	}
	return 0;
}