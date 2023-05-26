#include "SabreArray.h"
#include <fstream>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"


SabreArray::SabreArray() : 
	DetectorArray(), m_deadlayerEloss({14}, {28}, {1}, s_deadlayerThickness), 
	m_detectorEloss({14}, {28}, {1}, s_detectorThickness), m_degraderEloss({73}, {181}, {1}, s_degraderThickness)
{

	for(int i=0; i<s_nDets; i++)
		m_detectors.emplace_back(i, s_centerPhiList[i]*s_deg2rad, s_tilt*s_deg2rad, s_zOffset);
	//Only 0, 1, 4 valid in degrader land
	m_degradedDetectors[0] = true;
	m_degradedDetectors[1] = true;
	m_degradedDetectors[2] = false;
	m_degradedDetectors[3] = false;
	m_degradedDetectors[4] = true;
    //No degraded detectors
	// m_degradedDetectors[0] = false;
	// m_degradedDetectors[1] = false;
	// m_degradedDetectors[2] = false;
	// m_degradedDetectors[3] = false;
	// m_degradedDetectors[4] = false;

	//Choose who to look at right now. Usually switch on or off degraded/non-degraded.
	m_activeDetectors[0] = false;
	m_activeDetectors[1] = false;
	m_activeDetectors[2] = true;
	m_activeDetectors[3] = true;
	m_activeDetectors[4] = false;
}

SabreArray::~SabreArray() {}

void SabreArray::DrawDetectorSystem(const std::string& filename)
{
	std::ofstream output(filename);

	std::vector<double> ringxs, ringys, ringzs;
    std::vector<double> wedgexs, wedgeys, wedgezs;
	ROOT::Math::XYZPoint coords;
 	for(int i=0; i<s_nDets; i++)
	{
 		for(int j=0; j<m_detectors[i].GetNumberOfRings(); j++)
		{
 			for(int k=0; k<4; k++)
			{
 				coords = m_detectors[i].GetRingTiltCoords(j, k);
 				ringxs.push_back(coords.X());
 				ringys.push_back(coords.Y());
 				ringzs.push_back(coords.Z());
 			}
 		}
 	}

 	for(int i=0; i<s_nDets; i++)
	{
 		for(int j=0; j<m_detectors[i].GetNumberOfWedges(); j++)
		{
 			for(int k=0; k<4; k++)
			{
 				coords = m_detectors[i].GetWedgeTiltCoords(j, k);
 				wedgexs.push_back(coords.X());
 				wedgeys.push_back(coords.Y());
 				wedgezs.push_back(coords.Z());
 			}
 		}
 	}

	for(unsigned int i=0; i<ringxs.size(); i++)
		output<<ringxs[i]<<" "<<ringys[i]<<" "<<ringzs[i]<<std::endl;
	for(unsigned int i=0; i<wedgexs.size(); i++)
		output<<wedgexs[i]<<" "<<wedgeys[i]<<" "<<wedgezs[i]<<std::endl;

 	output.close();
}

double SabreArray::RunConsistencyCheck()
{
	double theta, phi;
	double npoints = 5.0*16.0*4.0;
	int count=0;
	ROOT::Math::XYZPoint corner;
 	for(auto& detector : m_detectors)
	{
 		for(int j=0; j<detector.GetNumberOfRings(); j++)
		{
 			for(int k=0; k<4; k ++) //Check corners
			{
				corner = detector.GetRingTiltCoords(j, k);
 				for(int i=0; i<5; i++)
				{
 					auto channels = m_detectors[i].GetTrajectoryRingWedge(corner.Theta(), corner.Phi());
 					if(channels.first != -1)
					{
 						count++;
						break;
 					}
 				}
 			}
 		}
 	}

 	return ((double)count)/npoints;
}

/*Returns if detected, as well as total energy deposited in SABRE*/
DetectorResult SabreArray::IsDetected(const Mask::Nucleus& nucleus)
{
	DetectorResult observation;
	if(nucleus.GetKE() <= s_energyThreshold)
		return observation;

	double thetaIncident;
	double ke = 0.0;
	for(auto& detector : m_detectors)
	{
		if(!m_activeDetectors[detector.GetDetectorID()])
			continue;

		auto channel = detector.GetTrajectoryRingWedge(nucleus.vec4.Theta(), nucleus.vec4.Phi());
		if(channel.first == -1 || channel.second == -1)
			continue;
		if(m_deadMap.IsDead(detector.GetDetectorID(), channel.first, 0) || m_deadMap.IsDead(detector.GetDetectorID(), channel.second, 1))
				break; //dead channel check

		observation.direction = detector.GetTrajectoryCoordinates(nucleus.vec4.Theta(), nucleus.vec4.Phi());
		thetaIncident = std::acos(observation.direction.Dot(detector.GetNormTilted())/(observation.direction.R()));

		//Energy loss
		ke = nucleus.GetKE();
		if(m_degradedDetectors[detector.GetDetectorID()])
			ke -= m_degraderEloss.GetEnergyLossTotal(nucleus.Z, nucleus.A, ke, M_PI - thetaIncident);
		ke -= m_deadlayerEloss.GetEnergyLossTotal(nucleus.Z, nucleus.A, ke, M_PI - thetaIncident);
		if(ke <= s_energyThreshold)
			break;

		observation.det_name = "SABRE"+std::to_string(detector.GetDetectorID());
		observation.energy_deposited = m_detectorEloss.GetEnergyLossTotal(nucleus.Z, nucleus.A, ke, M_PI - thetaIncident);
		observation.detectFlag = true;
		return observation;
	}

	observation.detectFlag = false;
	return observation;
}
