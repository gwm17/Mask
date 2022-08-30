#include "AnasenArray.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#include "TFile.h"
#include "TTree.h"

AnasenArray::AnasenArray() :
	DetectorArray(), m_detectorEloss({14}, {28}, {1}, s_detectorThickness)
{
	for(int i=0; i<s_nSX3PerBarrel; i++)
	{
		m_Ring1.emplace_back(s_barrelPhiList[i], s_barrel1Z, s_barrelRhoList[i]);
		m_Ring1[i].SetPixelSmearing(true);
		m_Ring2.emplace_back(s_barrelPhiList[i], s_barrel2Z, s_barrelRhoList[i]);
		m_Ring2[i].SetPixelSmearing(true);
	}
	for(int i=0; i<s_nQQQ; i++)
	{
		m_forwardQQQs.emplace_back(s_qqqPhiList[i], s_qqqZList[i]);
		m_forwardQQQs[i].SetSmearing(true);
		m_backwardQQQs.emplace_back(s_qqqPhiList[i], (-1.0)*s_qqqZList[i]);
		m_backwardQQQs[i].SetSmearing(true);
	}
}

AnasenArray::~AnasenArray() {}


void AnasenArray::DrawDetectorSystem(const std::string& filename)
{
	std::ofstream output(filename);

	std::vector<double> x, y, z;
	std::vector<double> cx, cy, cz;
	ROOT::Math::XYZPoint coords;
	for(int i=0; i<s_nSX3PerBarrel; i++)
	{
		for(int j=0; j<4; j++)
		{
			for(int k=0; k<4; k++)
			{
				coords = m_Ring1[i].GetRotatedFrontStripCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
				coords = m_Ring1[i].GetRotatedBackStripCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
			}
			coords = m_Ring1[i].GetHitCoordinates(j, 0);
			cx.push_back(coords.X());
			cy.push_back(coords.Y());
			cz.push_back(coords.Z());
		}
	}
	for(int i=0; i<s_nSX3PerBarrel; i++)
	{
		for(int j=0; j<4; j++)
		{
			for(int k=0; k<4; k++)
			{
				coords = m_Ring2[i].GetRotatedFrontStripCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
				coords = m_Ring2[i].GetRotatedBackStripCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
			}
			coords = m_Ring2[i].GetHitCoordinates(j, 0);
			cx.push_back(coords.X());
			cy.push_back(coords.Y());
			cz.push_back(coords.Z());
		}
	}
	for(int i=0; i<s_nQQQ; i++)
	{
		for(int j=0; j<16; j++)
		{
			for(int k=0; k<4; k++)
			{
				coords = m_forwardQQQs[i].GetRingCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
				coords = m_forwardQQQs[i].GetWedgeCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
			}
			for(int k=0; k<16; k++)
			{
				coords = m_forwardQQQs[i].GetHitCoordinates(j, k);
				cx.push_back(coords.X());
				cy.push_back(coords.Y());
				cz.push_back(coords.Z());
			}
		}
	}
	for(int i=0; i<s_nQQQ; i++)
	{
		for(int j=0; j<16; j++)
		{
			for(int k=0; k<4; k++)
			{
				coords = m_backwardQQQs[i].GetRingCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
				coords = m_backwardQQQs[i].GetWedgeCoordinates(j, k);
				x.push_back(coords.X());
				y.push_back(coords.Y());
				z.push_back(coords.Z());
			}
			for(int k=0; k<16; k++)
			{
				coords = m_backwardQQQs[i].GetHitCoordinates(j, k);
				cx.push_back(coords.X());
				cy.push_back(coords.Y());
				cz.push_back(coords.Z());
			}
		}
	}

	output<<"ANASEN Geometry File -- Coordinates for Detectors"<<std::endl;
	for(std::size_t i=0; i<x.size(); i++)
		output<<x[i]<<" "<<y[i]<<" "<<z[i]<<std::endl;
	for(std::size_t i=0; i<cx.size(); i++)
		output<<cx[i]<<" "<<cy[i]<<" "<<cz[i]<<std::endl;

	output.close();
}

double AnasenArray::RunConsistencyCheck()
{
	std::vector<ROOT::Math::XYZPoint> r1_points;
	std::vector<ROOT::Math::XYZPoint> r2_points;
	std::vector<ROOT::Math::XYZPoint> fqqq_points;
	std::vector<ROOT::Math::XYZPoint> bqqq_points;
	for(int i=0; i<s_nSX3PerBarrel; i++)
	{
		for(int j=0; j<4; j++)
			r1_points.push_back(m_Ring1[i].GetHitCoordinates(j, 0));
	}
	for(int i=0; i<s_nSX3PerBarrel; i++)
	{
		for(int j=0; j<4; j++)
			r2_points.push_back(m_Ring2[i].GetHitCoordinates(j, 0));
	}
	for(int i=0; i<s_nQQQ; i++)
	{
		for(int j=0; j<16; j++)
		{
			for(int k=0; k<16; k++)
				fqqq_points.push_back(m_forwardQQQs[i].GetHitCoordinates(j, k));
		}
	}
	for(int i=0; i<s_nQQQ; i++)
	{
		for(int j=0; j<16; j++)
		{
			for(int k=0; k<16; k++)
				bqqq_points.push_back(m_backwardQQQs[i].GetHitCoordinates(j, k));
		}
	}

	std::size_t npoints = r1_points.size() + r2_points.size() + fqqq_points.size() + bqqq_points.size();
	std::size_t count = 0;
	ROOT::Math::XYZPoint coords;
	for(auto& point : r1_points)
	{
		for(auto& sx3 : m_Ring1)
		{
			auto result = sx3.GetChannelRatio(point.Theta(), point.Phi());
			coords = sx3.GetHitCoordinates(result.front_strip_index, result.front_ratio);
			if(IsDoubleEqual(point.X(), coords.X()) && IsDoubleEqual(point.Y(), coords.Y()) && IsDoubleEqual(point.Z(), coords.Z()))
			{
				count++;
				break;
			}
		}
	}
	for(auto& point : r2_points)
	{
		for(auto& sx3 : m_Ring2)
		{
			auto result = sx3.GetChannelRatio(point.Theta(), point.Phi());
			coords = sx3.GetHitCoordinates(result.front_strip_index, result.front_ratio);
			if(IsDoubleEqual(point.X(), coords.X()) && IsDoubleEqual(point.Y(), coords.Y()) && IsDoubleEqual(point.Z(), coords.Z()))
			{
				count++;
				break;
			}
		}
	}
	for(auto& point : fqqq_points)
	{
		for(auto& qqq : m_forwardQQQs)
		{
			auto result = qqq.GetTrajectoryRingWedge(point.Theta(), point.Phi());
			coords = qqq.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.X(), coords.X()) && IsDoubleEqual(point.Y(), coords.Y()) && IsDoubleEqual(point.Z(), coords.Z()))
			{
				count++;
				break;
			}
		}
	}
	for(auto& point : bqqq_points)
	{
		for(auto& qqq : m_backwardQQQs)
		{
			auto result = qqq.GetTrajectoryRingWedge(point.Theta(), point.Phi());
			coords = qqq.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.X(), coords.X()) && IsDoubleEqual(point.Y(), coords.Y()) && IsDoubleEqual(point.Z(), coords.Z()))
			{
				count++;
				break;
			}
		}
	}

	double ratio = ((double)count)/((double)npoints);

	return ratio;

}

DetectorResult AnasenArray::IsRing1(Mask::Nucleus& nucleus)
{
	DetectorResult observation;
	double thetaIncident;
	for(int i=0; i<s_nSX3PerBarrel; i++)
	{
		auto result = m_Ring1[i].GetChannelRatio(nucleus.vec4.Theta(), nucleus.vec4.Phi());
		if(result.front_strip_index != -1 /*&& !dmap.IsDead(AnasenDetectorType::Barrel1, i, result.front_strip_index, AnasenDetectorSide::Front)*/
			&& !dmap.IsDead(AnasenDetectorType::Barrel1, i, result.back_strip_index, AnasenDetectorSide::Back)) 
		{
			observation.detectFlag = true;
			observation.direction = m_Ring1[i].GetHitCoordinates(result.front_strip_index, result.front_ratio);
			thetaIncident = std::acos(observation.direction.Dot(m_Ring1[i].GetNormRotated())/observation.direction.R());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			observation.energy_deposited = m_detectorEloss.GetEnergyLossTotal(nucleus.Z, nucleus.A, nucleus.GetKE(), thetaIncident);
			observation.det_name = "R1";
			return observation;
		}
	}

	return observation;
}

DetectorResult AnasenArray::IsRing2(Mask::Nucleus& nucleus)
{
	DetectorResult observation;
	double thetaIncident;
	for(int i=0; i<s_nSX3PerBarrel; i++)
	{
		auto result = m_Ring2[i].GetChannelRatio(nucleus.vec4.Theta(), nucleus.vec4.Phi());
		if(result.front_strip_index != -1 /*&& !dmap.IsDead(AnasenDetectorType::Barrel2, i, result.front_strip_index, AnasenDetectorSide::Front)*/
			&& !dmap.IsDead(AnasenDetectorType::Barrel2, i, result.back_strip_index, AnasenDetectorSide::Back)) 
		{
			observation.detectFlag = true;
			observation.direction = m_Ring2[i].GetHitCoordinates(result.front_strip_index, result.front_ratio);
			thetaIncident = std::acos(observation.direction.Dot(m_Ring2[i].GetNormRotated())/observation.direction.R());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			observation.energy_deposited = m_detectorEloss.GetEnergyLossTotal(nucleus.Z, nucleus.A, nucleus.GetKE(), thetaIncident);
			observation.det_name = "R2";
			return observation;
		}
	}

	return observation;
}

DetectorResult AnasenArray::IsQQQ(Mask::Nucleus& nucleus)
{
	DetectorResult observation;
	double thetaIncident;
	for(int i=0; i<s_nQQQ; i++)
	{
		auto result = m_forwardQQQs[i].GetTrajectoryRingWedge(nucleus.vec4.Theta(), nucleus.vec4.Phi());
		if(result.first != -1 /*&& !dmap.IsDead(AnasenDetectorType::FQQQ, i, result.first, AnasenDetectorSide::Front)*/ &&
			!dmap.IsDead(AnasenDetectorType::FQQQ, i, result.second, AnasenDetectorSide::Back)) 
		{
			observation.detectFlag = true;
			observation.direction = m_forwardQQQs[i].GetHitCoordinates(result.first, result.second);
			thetaIncident = std::acos(observation.direction.Dot(m_forwardQQQs[i].GetNorm())/observation.direction.R());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			observation.energy_deposited = m_detectorEloss.GetEnergyLossTotal(nucleus.Z, nucleus.A, nucleus.GetKE(), thetaIncident);
			observation.det_name = "FQQQ";
			return observation;
		}
	}

	
	for(int i=0; i<s_nQQQ; i++)
	{
		auto result = m_backwardQQQs[i].GetTrajectoryRingWedge(nucleus.vec4.Theta(), nucleus.vec4.Phi());
		if(result.first != -1 /*&& !dmap.IsDead(AnasenDetectorType::BQQQ, i, result.first, AnasenDetectorSide::Front)*/ &&
			!dmap.IsDead(AnasenDetectorType::BQQQ, i, result.second, AnasenDetectorSide::Back)) 
		{
			observation.detectFlag = true;
			observation.direction = m_backwardQQQs[i].GetHitCoordinates(result.first, result.second);
			thetaIncident = std::acos(observation.direction.Dot(m_backwardQQQs[i].GetNorm())/observation.direction.R());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			observation.energy_deposited = m_detectorEloss.GetEnergyLossTotal(nucleus.Z, nucleus.A, nucleus.GetKE(), thetaIncident);
			observation.det_name = "BQQQ";
			return observation;
		}
	}

	return observation;
}

DetectorResult AnasenArray::IsAnasen(Mask::Nucleus& nucleus)
{
	DetectorResult result;
	if(nucleus.GetKE() <= s_energyThreshold)
		return result;

	if(!result.detectFlag)
		result = IsRing1(nucleus);
	if(!result.detectFlag)
		result = IsRing2(nucleus);
	if(!result.detectFlag)
		result = IsQQQ(nucleus);
	return result;
}

void AnasenArray::CountCoincidences(const std::vector<Mask::Nucleus>& data, std::vector<int>& counts)
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

void AnasenArray::CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) {

	if(!dmap.IsValid())
	{
		std::cerr<<"Invalid Dead Channel Map at AnasenArray::CalculateEfficiency()! If you want to run with all possible";
		std::cerr<<"channels active, simply load an empty file."<<std::endl;
		return;
	}

	std::cout<<"----------ANASEN Efficiency Calculation----------"<<std::endl;
	std::cout<<"Loading in output from kinematics simulation: "<<inputname<<std::endl;
	std::cout<<"Running efficiency calculation..."<<std::endl;

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
			std::cerr<<"Bad reaction type at AnasenArray::CalculateEfficiency (given value: "<<counts.size()<<"). Quiting..."<<std::endl;
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
			DetectorResult result = IsAnasen(nucleus);
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
	for(unsigned int i=0; i<counts.size(); i++)
	{
		stats<<std::setw(10)<<i<<"|"<<std::setw(10)<<((double)counts[i]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	stats<<"Coincidence Efficiency"<<std::endl;
	stats<<"---------------------"<<std::endl;
	if(dataHandle->size() == 3)
	{
		stats<<std::setw(10)<<"1 + 2"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(dataHandle->size() == 4)
	{
		stats<<std::setw(10)<<"2 + 3"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/nentries)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(dataHandle->size() == 6)
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
	else if(dataHandle->size() == 8)
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