#ifndef ANASEN_EFFICIENCY_H
#define ANASEN_EFFICIENCY_H

#include <string>

#include "DetectorEfficiency.h"
#include "StripDetector.h"
#include "QQQDetector.h"
#include "Target.h"
#include "Nucleus.h"
#include "MaskFile.h"

struct DetectorResult {
	bool detectFlag = false;
	Mask::Vec3 direction;
	double energy_deposited = 0.0;
	std::string det_name = "";
};

class AnasenEfficiency : public DetectorEfficiency {
public:
	AnasenEfficiency();
	~AnasenEfficiency();
	void CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) override;
	void DrawDetectorSystem(const std::string& filename) override;
	double RunConsistencyCheck() override;

private:
	DetectorResult IsRing1(Mask::Nucleus& nucleus);
	DetectorResult IsRing2(Mask::Nucleus& nucleus);
	DetectorResult IsQQQ(Mask::Nucleus& nucleus);
	DetectorResult IsAnasen(Mask::Nucleus& nucleus);
	void CountCoincidences(const Mask::MaskFileData& data, std::vector<int>& counts, int rxn_type);

	std::vector<StripDetector> m_Ring1, m_Ring2;
	std::vector<QQQDetector> m_forwardQQQs;
	std::vector<QQQDetector> m_backwardQQQs;

	Mask::Target det_silicon;

	/**** ANASEN geometry constants *****/
	const int n_sx3_per_ring = 12;
	const int n_qqq = 4;
	const double sx3_length = 0.075;
	const double sx3_width = 0.04;
	const double barrel_gap = 0.0254;
	const double sx3_frame = 0.049; //0.049 is base gap due to frames
	const double ring1_z = sx3_length/2.0 + sx3_frame + barrel_gap/2.0;
	const double ring2_z = (-1.0)*(barrel_gap/2.0 + sx3_length/2.0);
	const double qqq_nom_z = 0.0125 + sx3_length + sx3_frame + barrel_gap/2.0;
	const double qqq_rinner = 0.0501;
	const double qqq_router = 0.0990;
	const double qqq_deltaphi = 1.52119;
	const double qqq_z[4] = {qqq_nom_z, qqq_nom_z, qqq_nom_z, qqq_nom_z};
	const double qqq_phi[4] = {5.49779, 0.785398, 2.35619, 3.92699};
	const double ring_rho[12] = {0.0890601, 0.0889871, 0.0890354, 0.0890247, 0.0890354, 0.0890354, 0.0890247, 0.0890354, 0.0890354, 0.0890247, 0.0890354, 0.0890354};
	const double ring_phi[12] = {0.785795, 0.262014, 6.02132, 5.49779, 4.97426, 4.45052, 3.92699, 3.40346, 2.87972, 2.35619, 1.83266, 1.30893};
	/*************************/

	static constexpr double threshold = 0.6; //MeV
	static constexpr double deg2rad = M_PI/180.0;
	static constexpr double si_thickness = 1000 * 1e-4 * 2.3926 * 1e6; //thickness in um -> eff thickness in ug/cm^2 for detector
};

#endif