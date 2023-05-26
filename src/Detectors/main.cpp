#include "DetectorApp.h"
#include "KinematicsExceptions.h"
#include <iostream>
#include <string>

int main(int argc, char** argv)
{

	if(argc != 2)
	{
		std::cerr<<"Incorrect number of commandline arguments! Returning."<<std::endl;
		return 1;
	}

	if(!Mask::EnforceDictionaryLinked())
	{
		std::cerr<<"This should be illegal!"<<std::endl;
		return 1;
	}

	try
	{
		DetectorApp app;
		if(!app.ParseConfig(argv[1]))
		{
			std::cerr << "Unable to load config file " << argv[1] << ". Shutting down." << std::endl;
			return 1;
		}

		app.Run();
	}
	catch(std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return 1;
	}
	return 0;
}
