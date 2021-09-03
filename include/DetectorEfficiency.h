#ifndef DETECTOREFFICIENCY_H
#define DETECTOREFFICIENCY_H

#include <string>
#include <cmath>

class DetectorEfficiency {
public:
	DetectorEfficiency() {};
	virtual ~DetectorEfficiency() {};

	virtual void CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) = 0;
	virtual void DrawDetectorSystem(const std::string& filename) = 0;
	virtual double RunConsistencyCheck() = 0;

protected:
	inline bool IsDoubleEqual(double x, double y) { return std::fabs(x-y) < epsilon ? true : false; };

	static constexpr double epsilon = 1.0e-6;
};

#endif