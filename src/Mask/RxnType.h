#ifndef RXNTYPE_H
#define RXNTYPE_H

#include <string>

namespace Mask {

	enum class RxnType
	{
		Decay,
		Reaction,
		None
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

	static RxnType StringToRxnType(const std::string& type)
	{
		if (type == "Decay")
			return RxnType::Decay;
		else if (type == "Reaction")
			return RxnType::Reaction;
		else
			return RxnType::None;
	}

	static std::string RxnTypeToString(RxnType type)
	{
		switch(type)
		{
			case RxnType::Decay: return "Decay";
			case RxnType::Reaction: return "Reaction";
			case RxnType::None: return "None";
			default: return "None";
		}
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

	static std::string RxnThetaTypeToString(RxnThetaType type)
	{
		switch(type)
		{
			case RxnThetaType::CenterOfMass: return "CenterOfMass";
			case RxnThetaType::Lab: return "Lab";
			case RxnThetaType::None: return "None";
			default: return "None";
		}
	}
}

#endif