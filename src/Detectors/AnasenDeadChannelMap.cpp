#include "AnasenDeadChannelMap.h"
#include <fstream>
#include <iostream>
#include <cmath>

AnasenDeadChannelMap::AnasenDeadChannelMap() :
	valid_flag(false)
{
	InitMap();
}

AnasenDeadChannelMap::AnasenDeadChannelMap(const std::string& filename) :
	valid_flag(false)
{
	InitMap();
	LoadMapfile(filename);
}

AnasenDeadChannelMap::~AnasenDeadChannelMap() {}

void AnasenDeadChannelMap::InitMap()
{
	for(int i=0; i<(bqqq_offset+4.0*nchannels_qqq); i++)
		dcMap[i] = false;
}

int AnasenDeadChannelMap::ConvertStringTypeIndexToOffset(const std::string& type, const std::string& index, const std::string& side)
{
	int channel_offset=0;
	if(type == "BARREL1B")
		channel_offset += 6*nchannels_sx3;
	else if (type == "BARREL2A")
		channel_offset += barrel2_offset;
	else if (type == "BARREL2B")
		channel_offset += barrel2_offset + 6*nchannels_sx3;
	else if (type == "FQQQ")
		channel_offset += fqqq_offset;
	else if (type == "BQQQ")
		channel_offset += bqqq_offset;

	if(index == "B")
		channel_offset += nchannels_sx3;
	else if(index == "C")
		channel_offset += 2.0*nchannels_sx3;
	else if(index == "D")
		channel_offset += 3.0*nchannels_sx3;
	else if(index == "E")
		channel_offset += 4.0*nchannels_sx3;
	else if(index == "F")
		channel_offset += 5.0*nchannels_sx3;
	else if(index == "1")
		channel_offset += nchannels_qqq;
	else if(index == "2")
		channel_offset += 2.0*nchannels_qqq;
	else if(index == "3")
		channel_offset += 3.0*nchannels_qqq;

	if(side == "BACK")
		channel_offset += nfronts_sx3;
	else if(side == "WEDGE")
		channel_offset += nfronts_qqq;

	return channel_offset;
}

void AnasenDeadChannelMap::LoadMapfile(const std::string& filename)
{
	valid_flag = false;
	std::ifstream input(filename);
	if(!input.is_open())
	{
		std::cerr<<"Unable to open input file at AnasenDeadChannelMap::LoadMapfile(). Exiting."<<std::endl;
		return;
	}

	std::string junk;
	std::string type, side, index;
	int channel;

	int dead_channel;
	while(input>>junk)
	{
		input>>type>>index>>side>>junk>>channel;

		if(side == "FRONT")
		{
			channel = std::floor(channel/2.0);
		}

		dead_channel = ConvertStringTypeIndexToOffset(type, index, side) + channel;
		
		dcMap[dead_channel] = true;
	}

	valid_flag = true;
}

const bool AnasenDeadChannelMap::IsDead(AnasenDetectorType type, int detIndex, int channel, AnasenDetectorSide side) const
{
	int channel_index=-1;
	switch(type)
	{
		case AnasenDetectorType::Barrel1:
		{
			channel_index = detIndex*nchannels_sx3 + side*nfronts_sx3 + channel;
			break;
		}
		case AnasenDetectorType::Barrel2:
		{
			channel_index = barrel2_offset + detIndex*nchannels_sx3 + side*nfronts_sx3 + channel;
			break;
		}
		case AnasenDetectorType::FQQQ:
		{
			channel_index = fqqq_offset + detIndex*nchannels_qqq + side*nfronts_qqq + channel;
			break;
		}
		case AnasenDetectorType::BQQQ:
		{
			channel_index = bqqq_offset + detIndex*nchannels_qqq + side*nfronts_qqq + channel;
			break;
		}
	}

	auto data = dcMap.find(channel_index);
	if(data == dcMap.end())
	{
		std::cerr<<"Bad channel index "<<channel_index<<" at AnasenDeadChannelMap::IsDead()"<<std::endl;
		return false;
	}

	return data->second;
}