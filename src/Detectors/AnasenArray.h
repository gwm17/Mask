#ifndef ANASEN_ARRAY_H
#define ANASEN_ARRAY_H

#include <string>

#include "DetectorArray.h"
#include "SX3Detector.h"
#include "QQQDetector.h"
#include "Mask/Target.h"
#include "Mask/Nucleus.h"
#include "AnasenDeadChannelMap.h"

class AnasenArray : public DetectorArray
{
public:
	AnasenArray();
	~AnasenArray();
	virtual DetectorResult IsDetected(const Mask::Nucleus& nucleus);
	virtual void DrawDetectorSystem(const std::string& filename) override;
	virtual double RunConsistencyCheck() override;
	virtual void SetDeadChannelMap(const std::string& filename) override { dmap.LoadMapfile(filename); }

private:
	DetectorResult IsRing1(const Mask::Nucleus& nucleus);
	DetectorResult IsRing2(const Mask::Nucleus& nucleus);
	DetectorResult IsQQQ(const Mask::Nucleus& nucleus);

	std::vector<SX3Detector> m_Ring1;
	std::vector<SX3Detector> m_Ring2;
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