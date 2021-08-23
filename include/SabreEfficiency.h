#ifndef SABREEFFICIENCY_H
#define SABREEFFICIENCY_H

#include "SabreDetector.h"
#include "Target.h"
#include "DeadChannelMap.h"
#include "Kinematics.h"
#include <THashTable.h>

class SabreEfficiency {
public:
	SabreEfficiency();
	~SabreEfficiency();
	inline void SetReactionType(int t) { m_rxn_type = t; };
    void SetDeadChannelMap(std::string& filename) { dmap.LoadMapfile(filename); };
	void CalculateEfficiency(const std::string& file);
    void DrawDetectorSystem(const std::string& filename);
    double RunConsistencyCheck();

private:
    void MyFill(THashTable* table, const std::string& name, const std::string& title, int bins, float min, float max, double val);
    void MyFill(THashTable* table, const std::string& name, const std::string& title, int binsx, float minx, float maxx, double valx, int binsy, float miny, float maxy, double valy);
	std::pair<bool,double> IsSabre(Mask::NucData* nucleus);
    void Run2Step(const std::string& filename);
	void Run3Step(const std::string& filename);
	void RunDecay(const std::string& filename);

	int m_rxn_type;
	std::vector<SabreDetector> detectors;
    
	Target deadlayer;
    Target sabre_eloss;
    DeadChannelMap dmap;


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
