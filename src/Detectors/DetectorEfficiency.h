#ifndef DETECTOREFFICIENCY_H
#define DETECTOREFFICIENCY_H

#include <string>
#include <cmath>

#include "Math/Point3D.h"

struct DetectorResult
{
	bool detectFlag = false;
	ROOT::Math::XYZPoint direction;
	double energy_deposited = 0.0;
	std::string det_name = "";
};

class DetectorEfficiency
{
public:
	DetectorEfficiency() {};
	virtual ~DetectorEfficiency() {};

	virtual void CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) = 0;
	virtual void DrawDetectorSystem(const std::string& filename) = 0;
	virtual double RunConsistencyCheck() = 0;

protected:
	inline bool IsDoubleEqual(double x, double y) { return std::fabs(x-y) < s_epsilon ? true : false; };

	static constexpr double s_epsilon = 1.0e-6;
};

#endif