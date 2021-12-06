#ifndef ANASENDEADCHANNELMAP_H
#define ANASENDEADCHANNELMAP_H

#include <unordered_map>
#include <string>

enum class AnasenDetectorType
{
	Barrel1,
	Barrel2,
	FQQQ,
	BQQQ
};

enum AnasenDetectorSide
{
	Front=0,
	Back=1
};

class AnasenDeadChannelMap
{
public:
	AnasenDeadChannelMap();
	AnasenDeadChannelMap(const std::string& filename);
	~AnasenDeadChannelMap();
	void LoadMapfile(const std::string& filename);
	inline const bool IsValid() const { return valid_flag; }
	const bool IsDead(AnasenDetectorType type, int detIndex, int channel, AnasenDetectorSide side) const;
private:
	void InitMap();
	int ConvertStringTypeIndexToOffset(const std::string& type, const std::string& index, const std::string& side);
	bool valid_flag;
	std::unordered_map<int, bool> dcMap;

	/*
		global channel calculated as 
		detTypeOffset + detIndex*total + detSideOffset + channel
	*/

	const int nfronts_sx3 = 4;
	const int nfronts_qqq = 16;
	const int nchannels_sx3 = 8;
	const int nchannels_qqq = 32;

	const int barrel1_offset = 0;
	const int barrel2_offset = 12*nchannels_sx3;
	const int fqqq_offset = barrel2_offset + 12*nchannels_sx3;
	const int bqqq_offset = fqqq_offset + 4*nchannels_qqq;

};

#endif