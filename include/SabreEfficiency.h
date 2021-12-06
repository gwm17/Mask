#ifndef SABREEFFICIENCY_H
#define SABREEFFICIENCY_H

#include "DetectorEfficiency.h"
#include "SabreDetector.h"
#include "Target.h"
#include "SabreDeadChannelMap.h"
#include "Nucleus.h"

class SabreEfficiency : public DetectorEfficiency {
public:
	SabreEfficiency();
	~SabreEfficiency();
    void SetDeadChannelMap(std::string& filename) { dmap.LoadMapfile(filename); };
	void CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) override;
    void DrawDetectorSystem(const std::string& filename) override;
    double RunConsistencyCheck() override;

private:
	std::pair<bool,double> IsSabre(Mask::Nucleus& nucleus);

	std::vector<SabreDetector> detectors;
    
	Mask::Target deadlayer;
    Mask::Target sabre_eloss;
    SabreDeadChannelMap dmap;

	//Sabre constants
	const double INNER_R = 0.0326;
    const double OUTER_R = 0.1351;
    const double TILT = 40.0;
    const double DIST_2_TARG = -0.1245;
    const double PHI_COVERAGE = 54.4; //delta phi for each det
    const double PHI0 = 234.0; //center phi values for each det in array
    const double PHI1 = 162.0; //# is equal to detID in channel map
    const double PHI2 = 306.0;
    const double PHI3 = 18.0;
    const double PHI4 = 90.0;
    const double DEG2RAD = M_PI/180.0;
    static constexpr double DEADLAYER_THIN = 50 * 1e-7 * 2.3296 * 1e6; // ug/cm^2 (50 nm thick * density)
    static constexpr double SABRE_THICKNESS = 500 * 1e-4 * 2.3926 * 1e6; // ug/cm^2 (500 um thick * density)

    const double ENERGY_THRESHOLD = 0.2; //in MeV

};

#endif
