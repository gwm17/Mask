#ifndef ANASEN_EFFICIENCY_H
#define ANASEN_EFFICIENCY_H

#include <string>

#include "DetectorEfficiency.h"
#include "StripDetector.h"
#include "QQQDetector.h"

class AnasenEfficiency : public DetectorEfficiency {
public:
	AnasenEfficiency();
	~AnasenEfficiency();
	void CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) override;
	void DrawDetectorSystem(const std::string& filename) override;
	double RunConsistencyCheck() override;

private:
	bool IsRing1(double theta, double phi);
	bool IsRing2(double theta, double phi);
	bool IsQQQ(double theta, double phi);

	std::vector<StripDetector> m_Ring1, m_Ring2;
	std::vector<QQQDetector> m_forwardQQQs;
	std::vector<QQQDetector> m_backwardQQQs;

	/**** ANASEN geometry constants *****/
	const int n_sx3_per_ring = 12;
	const int n_qqq = 4;
	const double sx3_length = 0.075;
	const double sx3_width = 0.04;
	const double barrel_gap = 0.013 + 0.049; //0.049 is base gap due to frames
	const double ring1_z = sx3_length/2.0 + barrel_gap/2.0;
	//const double ring2_z = -0.124 + sx3_length/2.0 + 0.0245 - barrel_gap/2.0;
	const double qqq_nom_z = 0.025 + sx3_length + 0.0245 + barrel_gap/2.0;
	const double qqq_rinner = 0.0501;
	const double qqq_router = 0.0990;
	const double qqq_deltaphi = 1.52119;
	const double qqq_z[4] = {qqq_nom_z, qqq_nom_z - 0.00828, qqq_nom_z, qqq_nom_z};
	const double qqq_phi[4] = {5.49779, 0.785398, 2.35619, 3.92699};
	const double ring_rho[12] = {0.0890601, 0.0889871, 0.0890354, 0.0890247, 0.0890354, 0.0890354, 0.0890247, 0.0890354, 0.0890354, 0.0890247, 0.0890354, 0.0890354};
	const double ring_phi[12] = {0.785795, 0.262014, 6.02132, 5.49779, 4.97426, 4.45052, 3.92699, 3.40346, 2.87972, 2.35619, 1.83266, 1.30893};
	/*************************/

	static constexpr double threshold = 0.2; //MeV
	static constexpr double deg2rad = M_PI/180.0;
};

#endif