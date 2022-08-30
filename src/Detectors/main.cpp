#include "SabreArray.h"
#include "AnasenArray.h"
#include "KinematicsExceptions.h"
#include <iostream>
#include <string>

int main(int argc, char** argv)
{

	if(argc != 4)
	{
		std::cerr<<"Incorrect number of commandline arguments! Returning."<<std::endl;
		return 1;
	}

	if(!Mask::EnforceDictionaryLinked())
	{
		std::cerr<<"This should be illegal!"<<std::endl;
		return 1;
	}

	std::string inputname = argv[1];
	std::string outputname = argv[2];
	std::string statsname = argv[3];

	
	SabreArray sabre;
	std::string mapfile = "./etc/sabreDeadChannels_May2022.txt";
	sabre.SetDeadChannelMap(mapfile);
	sabre.CalculateEfficiency(inputname, outputname, statsname);
	//std::cout<<"Running consistency check(1=success): "<<sabre.RunConsistencyCheck()<<std::endl;
	//sabre.DrawDetectorSystem("/data1/gwm17/10B3He/Feb2021/simulation/SABREGeo.txt");
	

	/*
 	AnasenArray anasen;
 	std::string mapfile = "./etc/AnasenDeadChannels.txt";
 	anasen.SetDeadChannelMap(mapfile);
 	anasen.CalculateEfficiency(inputname, outputname, statsname);
 	//std::cout<<"Running consistency check(1=success): "<<anasen.RunConsistencyCheck()<<std::endl;
 	//anasen.DrawDetectorSystem("/data1/gwm17/TRIUMF_7Bed/simulation/ANASENGeo_centered_target_targetGap_BackQQQ_fixedZ.txt");
	*/
	return 0;
}
