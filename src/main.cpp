#include <iostream>
#include <string>
#include "Stopwatch.h"
#include "MaskFile.h"
#include "Kinematics.h"
#include "KinematicsExceptions.h"

int main(int argc, char** argv) {
	if(argc<2) {
		std::cerr<<"Incorrect number of arguments!"<<std::endl;
		return 1;
	}
	Stopwatch sw;


	Mask::Kinematics calculator;
	sw.Start();
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
	sw.Stop();
	std::cout<<"Time elapsed(seconds): "<<sw.GetElapsedSeconds()<<std::endl;
	
	Mask::MaskFile input(calculator.GetOutputName(), Mask::MaskFile::FileType::read);
	Mask::MaskFileData data;
	Mask::MaskFileHeader header = input.ReadHeader();

	std::cout<<"Header Found -- rxn type: "<<header.rxn_type<<" nsamples: "<<header.nsamples;
	std::cout<<std::endl;

	int counter=0;
	while(!data.eof) {
		data = input.ReadData();
		for(unsigned int i=0; i<data.E.size(); i++)
			counter++;
			//std::cout<<"Data Found -- E: "<<data.E[i]<<" KE: "<<data.KE[i]<<" p: "<<data.p[i]<<" theta: "<<data.theta[i]<<" phi: "<<data.phi[i]<<std::endl;
	}
	std::cout<<"events found: "<<counter<<std::endl;
	input.Close();
	std::cout<<"File closed."<<std::endl;

	return 0;

}