#include "DeadChannelMap.h"
#include <fstream>

DeadChannelMap::DeadChannelMap() :
	validFlag(false)
{
	int maxchan = 5*24+1*16+7;
	for(int i=0; i<maxchan; i++) {
		dcMap[i] = false;
	}
}

DeadChannelMap::DeadChannelMap(std::string& name) :
	validFlag(false)
{
	LoadMapfile(name);
}

DeadChannelMap::~DeadChannelMap() {}

void DeadChannelMap::LoadMapfile(std::string& name) {
	std::ifstream input(name);
	if(!input.is_open()) {
		std::cerr<<"Unable to load dead channels, file: "<<name<<" could not be opened"<<std::endl;
		validFlag = false;
		return;
	}

	

	std::string junk, rw;
	int detID, channel;
	int this_channel;

	std::getline(input, junk);
	while(input>>detID) {
		input>>rw>>channel;
		if(rw == "RING") {
			this_channel = detID*24+RING*16+channel;
			dcMap[this_channel] = true;
		} else if(rw == "WEDGE") {
			this_channel = detID*24+WEDGE*16+channel;
			dcMap[this_channel] = true;
		} else {
			std::cerr<<"Invalid ring/wedge designation at DeadChannelMap"<<std::endl;
		}
	}

	validFlag = true;
}

bool DeadChannelMap::IsDead(int detID, int channel, int ringwedgeFlag) {
	if(!IsValid() || (ringwedgeFlag != 0 && ringwedgeFlag != 1)) {
		std::cerr<<"Error at DeadChannelMap IsDead(), bad input parameters"<<std::endl;
		return false;
	}

	int this_channel = detID*24+ringwedgeFlag*16+channel;

	auto iter = dcMap.find(this_channel);
	if(iter != dcMap.end()) {
		return iter->second;
	} else {
		std::cerr<<"Error at DeadChannelMap IsDead(), invalid channel: "<<this_channel<<std::endl;
		std::cerr<<"detID: "<<detID<<" ringwedgeFlag: "<<ringwedgeFlag<<" channel: "<<channel<<std::endl;
		return false;
	}
}