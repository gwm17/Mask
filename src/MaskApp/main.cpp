#include <iostream>
#include <string>
#include "Stopwatch.h"
#include "MaskFile.h"
#include "MaskApp.h"
#include "KinematicsExceptions.h"

int main(int argc, char** argv) {
	if(argc<2) {
		std::cerr<<"Incorrect number of arguments!"<<std::endl;
		return 1;
	}
	Stopwatch sw;


	Mask::MaskApp calculator;
	sw.Start();
	try {
		if(!calculator.LoadConfig(argv[1])) {
			std::cerr<<"Unable to read input file!"<<std::endl;
			return 1;
		}
		calculator.Run();
	} catch(const std::exception& e) {
		std::cerr<<"Exception caught! Information: "<<e.what()<<std::endl;
		std::cerr<<"Terminating process."<<std::endl;
		return 1;
	}
	sw.Stop();
	std::cout<<"Time elapsed(seconds): "<<sw.GetElapsedSeconds()<<std::endl;
	
	return 0;

}