#ifndef SABREDEADCHANNELMAP_H
#define SABREDEADCHANNELMAP_H

#include <string>
#include <iostream>
#include <unordered_map>

class SabreDeadChannelMap 
{
public:
	SabreDeadChannelMap();
	SabreDeadChannelMap(const std::string& name);
	~SabreDeadChannelMap();
	void LoadMapfile(const std::string& name);
	bool IsDead(int det, int channel, int ringwedgeFlag);
	bool IsValid() { return validFlag; };

private:
	std::unordered_map<int, bool> dcMap;
	bool validFlag;
	static constexpr int RING = 0;
	static constexpr int WEDGE = 1;

	/*
	Channel identifier calculated like detector*24 + 16*(RING or WEDGE) + channel
	Example:
	Detector 1 ring 15
	1*24 + 16*0 + 15 = 39
	Detector 1 wedge 0
	1*24 + 16*1 +0= 40
	*/
};

#endif