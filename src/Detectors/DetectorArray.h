#ifndef DETECTOR_ARRAY_H
#define DETECTOR_ARRAY_H

#include <string>
#include <cmath>

#include "Math/Point3D.h"
#include "Mask/Nucleus.h"

struct DetectorResult
{
	bool detectFlag = false;
	ROOT::Math::XYZPoint direction;
	double energy_deposited = 0.0;
	std::string det_name = "";
};

enum class ArrayType
{
	None,
	Anasen,
	Sabre
};

class DetectorArray
{
public:
	DetectorArray() {};
	virtual ~DetectorArray() {};

	virtual DetectorResult IsDetected(const Mask::Nucleus& nucleus) = 0;
	virtual void DrawDetectorSystem(const std::string& filename) = 0;
	virtual double RunConsistencyCheck() = 0;
	virtual void SetDeadChannelMap(const std::string& filename) = 0;

protected:
	bool IsDoubleEqual(double x, double y) { return std::fabs(x-y) < s_epsilon ? true : false; };

	static constexpr double s_epsilon = 1.0e-6;
};

DetectorArray* CreateDetectorArray(ArrayType type);

std::string ArrayTypeToString(ArrayType type);

ArrayType StringToArrayType(const std::string& value);

#endif