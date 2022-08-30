#ifndef ANASEN_EFFICIENCY_H
#define ANASEN_EFFICIENCY_H

#include <string>

#include "DetectorEfficiency.h"
#include "SX3Detector.h"
#include "QQQDetector.h"
#include "Target.h"
#include "Nucleus.h"
#include "AnasenDeadChannelMap.h"

class AnasenEfficiency : public DetectorEfficiency
{
public:
	AnasenEfficiency();
	~AnasenEfficiency();
	void CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) override;
	void DrawDetectorSystem(const std::string& filename) override;
	double RunConsistencyCheck() override;
	inline void SetDeadChannelMap(const std::string& filename) { dmap.LoadMapfile(filename); }

private:
	DetectorResult IsRing1(Mask::Nucleus& nucleus);
	DetectorResult IsRing2(Mask::Nucleus& nucleus);
	DetectorResult IsQQQ(Mask::Nucleus& nucleus);
	DetectorResult IsAnasen(Mask::Nucleus& nucleus);
	void CountCoincidences(const std::vector<Mask::Nucleus>& data, std::vector<int>& counts);

	std::vector<StripDetector> m_Ring1, m_Ring2;
	std::vector<QQQDetector> m_forwardQQQs;
	std::vector<QQQDetector> m_backwardQQQs;

	Mask::Target m_detectorEloss;

	AnasenDeadChannelMap dmap;

	/**** ANASEN geometry constants *****/
	static constexpr int s_nSX3PerBarrel = 12;
	static constexpr int s_nQQQ = 4;
	static constexpr double s_sx3Length = 0.075;
	static constexpr double s_barrelGap = 0.0254;
	static constexpr double s_sx3FrameGap = 0.049; //0.049 is base gap due to frames
	static constexpr double s_barrel1Z = s_sx3Length/2.0 + s_sx3FrameGap + s_barrelGap/2.0;
	static constexpr double s_barrel2Z = (-1.0)*(s_barrelGap/2.0 + s_sx3Length/2.0);
	static constexpr double s_qqqZ = 0.0125 + s_sx3Length + s_sx3FrameGap + s_barrelGap/2.0;
	static constexpr double s_qqqZList[4] = {s_qqqZ, s_qqqZ, s_qqqZ, s_qqqZ};
	static constexpr double s_qqqPhiList[4] = {5.49779, 0.785398, 2.35619, 3.92699};
	static constexpr double s_barrelRhoList[12] = {0.0890601, 0.0889871, 0.0890354, 0.0890247, 0.0890354, 0.0890354, 0.0890247,
												 0.0890354, 0.0890354, 0.0890247, 0.0890354, 0.0890354};
	static constexpr double s_barrelPhiList[12] = {4.97426, 5.49739, 6.02132, 0.261868, 0.785398, 1.30893, 1.83266, 2.35619, 2.87972, 
												   3.40346, 3.92699, 4.45052};
	/*************************/

	static constexpr double s_energyThreshold = 0.6; //MeV
	static constexpr double s_deg2rad = M_PI/180.0;
	static constexpr double s_detectorThickness = 1000 * 1e-4 * 2.3926 * 1e6; //thickness in um -> eff thickness in ug/cm^2 for detector
};

#endif