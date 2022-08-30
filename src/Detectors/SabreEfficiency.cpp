#include "SabreEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"


SabreEfficiency::SabreEfficiency() : 
	DetectorEfficiency(), m_deadlayerEloss({14}, {28}, {1}, s_deadlayerThickness), 
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

	//Choose who to look at right now. Usually switch on or off degraded/non-degraded.
	m_activeDetectors[0] = false;
	m_activeDetectors[1] = false;
	m_activeDetectors[2] = true;
	m_activeDetectors[3] = true;
	m_activeDetectors[4] = false;
}

SabreEfficiency::~SabreEfficiency() {}

void SabreEfficiency::CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname)
{
	std::cout<<"----------SABRE Efficiency Calculation----------"<<std::endl;
	std::cout<<"Loading in output from kinematics simulation: "<<inputname<<std::endl;
	std::cout<<"Running efficiency calculation..."<<std::endl;

	if(!m_deadMap.IsValid())
	{
		std::cerr<<"Unable to run SABRE Efficiency without a dead channel map."<<std::endl;
		std::cerr<<"If you have no dead channels, simply make a file that's empty"<<std::endl;
		std::cerr<<"Exiting."<<std::endl;
		std::cout<<"---------------------------------------------"<<std::endl;
	}


	TFile* input = TFile::Open(inputname.c_str(), "READ");
	TFile* output = TFile::Open(outputname.c_str(), "RECREATE");
	std::ofstream stats(statsname);
	stats<<std::setprecision(5);

	TTree* intree = (TTree*) input->Get("SimTree");
	std::vector<Mask::Nucleus>* dataHandle = new std::vector<Mask::Nucleus>();
	intree->SetBranchAddress("nuclei", &dataHandle);

	output->cd();
	TTree* outtree = new TTree("SimTree", "SimTree");
	outtree->Branch("nuclei", dataHandle);

	input->cd();

	stats<<"Efficiency statistics for data from "<<inputname<<" using the ANASEN geometry"<<std::endl;
	stats<<"Given in order of target=0, projectile=1, ejectile=2, residual=3, .... etc."<<std::endl;

	intree->GetEntry(1);
	std::vector<int> counts;
	std::vector<int> coinc_counts;
	counts.resize(dataHandle->size());
	switch(counts.size())
	{
		case 3: coinc_counts.resize(1, 0); break;
		case 4: coinc_counts.resize(1, 0); break;
		case 6: coinc_counts.resize(4, 0); break;
		case 8: coinc_counts.resize(11, 0); break;
		default:
		{
			std::cerr<<"Bad reaction type at AnasenEfficiency::CalculateEfficiency (given value: "<<counts.size()<<"). Quiting..."<<std::endl;
			input->Close();
			output->Close();
			stats.close();
			return;
		}
	}

	uint64_t nentries = intree->GetEntries();
	uint64_t percent5 = nentries*0.05;
	uint64_t count = 0;
	uint64_t npercent = 0;

	Mask::Nucleus nucleus;
	for(uint64_t i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		if(++count == percent5)
		{//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		for(std::size_t j=0; j<dataHandle->size(); j++)
		{
			Mask::Nucleus& nucleus = (*dataHandle)[j];
			DetectorResult result = IsSabre(nucleus);
			if(result.detectFlag)
			{
				nucleus.isDetected = true;
				nucleus.detectedKE = result.energy_deposited;
				nucleus.detectedTheta = result.direction.Theta();
				nucleus.detectedPhi = result.direction.Phi();
				counts[j]++;
			}
		}

		CountCoincidences(*dataHandle, coinc_counts);

		outtree->Fill();
	}
	input->Close();
	output->cd();
	outtree->Write(outtree->GetName(), TObject::kOverwrite);
	output->Close();

	delete dataHandle;

	stats<<std::setw(10)<<"Index"<<"|"<<std::setw(10)<<"Efficiency"<<std::endl;
	stats<<"---------------------"<<std::endl;
	for(unsigned int i=0; i<counts.size(); i++) {
		stats<<std::setw(10)<<i<<"|"<<std::setw(10)<<((double)counts[i]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	stats<<"Coincidence Efficiency"<<std::endl;
	stats<<"---------------------"<<std::endl;
	if(counts.size() == 3)
	{
		stats<<std::setw(10)<<"1 + 2"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(counts.size() == 4)
	{
		stats<<std::setw(10)<<"2 + 3"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(counts.size() == 6)
	{
		stats<<std::setw(10)<<"2 + 4"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 5"<<"|"<<std::setw(10)<<((double)coinc_counts[1]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 5"<<"|"<<std::setw(10)<<((double)coinc_counts[2]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 5"<<"|"<<std::setw(10)<<((double)coinc_counts[3]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(counts.size() == 8)
	{
		stats<<std::setw(10)<<"2 + 4"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 6"<<"|"<<std::setw(10)<<((double)coinc_counts[1]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[2]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 6"<<"|"<<std::setw(10)<<((double)coinc_counts[3]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[4]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[5]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 6"<<"|"<<std::setw(10)<<((double)coinc_counts[6]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[7]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[8]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[9]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[10]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	stats.close();

	std::cout<<std::endl;
	std::cout<<"Complete."<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
}

void SabreEfficiency::DrawDetectorSystem(const std::string& filename)
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

 	output<<"SABRE Geometry File -- Coordinates for Detectors"<<std::endl;
	output<<"Edges: x y z"<<std::endl;
	for(unsigned int i=0; i<ringxs.size(); i++)
		output<<ringxs[i]<<" "<<ringys[i]<<" "<<ringzs[i]<<std::endl;
	for(unsigned int i=0; i<wedgexs.size(); i++)
		output<<wedgexs[i]<<" "<<wedgeys[i]<<" "<<wedgezs[i]<<std::endl;

 	output.close();
}

double SabreEfficiency::RunConsistencyCheck()
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
 					}
 				}
 			}
 		}
 	}

 	return ((double)count)/npoints;
}

/*Returns if detected, as well as total energy deposited in SABRE*/
DetectorResult SabreEfficiency::IsSabre(Mask::Nucleus& nucleus)
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

		observation.detectFlag = true;
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
		return observation;
	}

	observation.detectFlag = false;
	return observation;
}

void SabreEfficiency::CountCoincidences(const std::vector<Mask::Nucleus>& data, std::vector<int>& counts)
{
	if (data.size() == 3 && data[1].isDetected && data[2].isDetected)
	{
		counts[0]++;
	}
	else if (data.size() == 4 && data[2].isDetected && data[3].isDetected)
	{
		counts[0]++;
	}
	else if(data.size() == 6)
	{
		if(data[2].isDetected && data[4].isDetected) 
		{
			counts[0]++;
		}
		if(data[2].isDetected && data[5].isDetected)
		{
			counts[1]++;
		}
		if(data[4].isDetected && data[5].isDetected)
		{
			counts[2]++;
		}
		if(data[2].isDetected && data[4].isDetected && data[5].isDetected)
		{
			counts[3]++;
		}
	}
	else if(data.size() == 8)
	{
		if(data[2].isDetected && data[4].isDetected) 
		{
			counts[0]++;
		}
		if(data[2].isDetected && data[6].isDetected)
		{
			counts[1]++;
		}
		if(data[2].isDetected && data[7].isDetected)
		{
			counts[2]++;
		}
		if(data[4].isDetected && data[6].isDetected)
		{
			counts[3]++;
		}
		if(data[4].isDetected && data[7].isDetected)
		{
			counts[4]++;
		}
		if(data[6].isDetected && data[7].isDetected)
		{
			counts[5]++;
		}
		if(data[2].isDetected && data[4].isDetected && data[6].isDetected)
		{
			counts[6]++;
		}
		if(data[2].isDetected && data[4].isDetected && data[7].isDetected)
		{
			counts[7]++;
		}
		if(data[2].isDetected && data[6].isDetected && data[7].isDetected)
		{
			counts[8]++;
		}
		if(data[4].isDetected && data[6].isDetected && data[7].isDetected)
		{
			counts[9]++;
		}
		if(data[2].isDetected && data[4].isDetected && data[6].isDetected && data[7].isDetected)
		{
			counts[10]++;
		}
	}
}
