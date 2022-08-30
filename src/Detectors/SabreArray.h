#ifndef SABRE_ARRAY_H
#define SABRE_ARRAY_H

#include "DetectorArray.h"
#include "SabreDetector.h"
#include "Target.h"
#include "SabreDeadChannelMap.h"
#include "Mask/Nucleus.h"

class SabreArray : public DetectorArray
{
public:
	SabreArray();
	~SabreArray();
    void SetDeadChannelMap(const std::string& filename) { m_deadMap.LoadMapfile(filename); };
	void CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) override;
    void DrawDetectorSystem(const std::string& filename) override;
    double RunConsistencyCheck() override;

private:
	DetectorResult IsSabre(Mask::Nucleus& nucleus);
    void CountCoincidences(const std::vector<Mask::Nucleus>& data, std::vector<int>& counts);

	std::vector<SabreDetector> m_detectors;
    
	Mask::Target m_deadlayerEloss;
    Mask::Target m_detectorEloss;
    Mask::Target m_degraderEloss;
    SabreDeadChannelMap m_deadMap;

    bool m_activeDetectors[5];
    bool m_degradedDetectors[5];

	//Sabre constants
    static constexpr double s_tilt = 40.0;
    static constexpr double s_zOffset = -0.1245;
    static constexpr int s_nDets = 5;
    static constexpr double s_centerPhiList[s_nDets] = {306.0, 18.0, 234.0, 162.0, 90.0};
    static constexpr double s_deg2rad = M_PI/180.0;
    static constexpr double s_deadlayerThickness = 50 * 1e-7 * 2.3296 * 1e6; // ug/cm^2 (50 nm thick * density)
    static constexpr double s_detectorThickness = 500 * 1e-4 * 2.3926 * 1e6; // ug/cm^2 (500 um thick * density)
    static constexpr double s_degraderThickness = 70.0 * 1.0e-4 * 16.69 * 1e6; //tantalum degrader (70 um thick) 

    static constexpr double s_energyThreshold = 0.25; //in MeV

};

#endif
