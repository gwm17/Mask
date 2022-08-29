#ifndef RXNTYPE_H
#define RXNTYPE_H

#include <string>

namespace Mask {

	enum class RxnType
	{
		PureDecay=0,
		OneStepRxn=1,
		TwoStepRxn=2,
		ThreeStepRxn=3,
		None=4
	};

	enum class RxnThetaType
	{
		CenterOfMass,
		Lab,
		None
	};

	enum RxnSize
	{
		DecaySize = 2,
		OneStepSize = 3,
		TwoStepSize = 4,
		ThreeStepSize = 5
	};

	static RxnType GetRxnTypeFromString(const std::string& type_str)
	{
		if (type_str == "PureDecay")
			return RxnType::PureDecay;
		else if (type_str == "OneStepRxn")
			return RxnType::OneStepRxn;
		else if (type_str == "TwoStepRxn")
			return RxnType::TwoStepRxn;
		else if (type_str == "ThreeStepRxn")
			return RxnType::ThreeStepRxn;
		else
			return RxnType::None;
	}

	static std::string GetStringFromRxnType(RxnType type)
	{
		switch(type)
		{
			case RxnType::PureDecay: return "PureDecay";
			case RxnType::OneStepRxn: return "OneStepRxn";
			case RxnType::TwoStepRxn: return "TwoStepRxn";
			case RxnType::ThreeStepRxn: return "ThreeStepRxn";
			case RxnType::None : return "None";
		}

		return "None";
	}

	static RxnType GetRxnTypeFromInt(uint32_t value)
	{
		return static_cast<RxnType>(value);
	}

	static uint32_t GetIntFromRxnType(RxnType type)
	{
		return static_cast<uint32_t>(type);
	}

	static RxnThetaType StringToRxnThetaType(const std::string& type)
	{
		if(type == "CenterOfMass")
		{
			return RxnThetaType::CenterOfMass;
		}
		else if(type == "Lab")
		{
			return RxnThetaType::Lab;
		}
		else
		{
			return RxnThetaType::None;
		}
	}
}

#endif