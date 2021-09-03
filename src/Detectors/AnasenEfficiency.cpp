#include "AnasenEfficiency.h"
#include "Kinematics.h"
#include "MaskFile.h"
#include <fstream>
#include <iomanip>

AnasenEfficiency::AnasenEfficiency() :
	DetectorEfficiency()
{
	for(int i=0; i<n_sx3_per_ring; i++) {
		m_Ring1.emplace_back(4, sx3_length, sx3_width, ring_phi[i], ring1_z, ring_rho[i]);
		m_Ring2.emplace_back(4, sx3_length, sx3_width, ring_phi[i], -1.0*ring1_z, ring_rho[i]);
	}
	for(int i=0; i<n_qqq; i++) {
		m_forwardQQQs.emplace_back(qqq_rinner, qqq_router, qqq_deltaphi, qqq_phi[i], qqq_z[i]);
		m_backwardQQQs.emplace_back(qqq_rinner, qqq_router, qqq_deltaphi, qqq_phi[i], (-1.0)*qqq_z[i]);
	}
}

AnasenEfficiency::~AnasenEfficiency() {}


void AnasenEfficiency::DrawDetectorSystem(const std::string& filename) {
	std::ofstream output(filename);

	std::vector<double> x, y, z;
	std::vector<double> cx, cy, cz;
	Mask::Vec3 coords;
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			for(int k=0; k<4; k++) {
				coords = m_Ring1[i].GetRotatedFrontStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_Ring1[i].GetRotatedBackStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			coords = m_Ring1[i].GetHitCoordinates(j, 0);
			cx.push_back(coords.GetX());
			cy.push_back(coords.GetY());
			cz.push_back(coords.GetZ());
		}
	}
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			for(int k=0; k<4; k++) {
				coords = m_Ring2[i].GetRotatedFrontStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_Ring2[i].GetRotatedBackStripCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			coords = m_Ring2[i].GetHitCoordinates(j, 0);
			cx.push_back(coords.GetX());
			cy.push_back(coords.GetY());
			cz.push_back(coords.GetZ());
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<4; k++) {
				coords = m_forwardQQQs[i].GetRingCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_forwardQQQs[i].GetWedgeCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			for(int k=0; k<16; k++) {
				coords = m_forwardQQQs[i].GetHitCoordinates(j, k);
				cx.push_back(coords.GetX());
				cy.push_back(coords.GetY());
				cz.push_back(coords.GetZ());
			}
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<4; k++) {
				coords = m_backwardQQQs[i].GetRingCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
				coords = m_backwardQQQs[i].GetWedgeCoordinates(j, k);
				x.push_back(coords.GetX());
				y.push_back(coords.GetY());
				z.push_back(coords.GetZ());
			}
			for(int k=0; k<16; k++) {
				coords = m_backwardQQQs[i].GetHitCoordinates(j, k);
				cx.push_back(coords.GetX());
				cy.push_back(coords.GetY());
				cz.push_back(coords.GetZ());
			}
		}
	}

	output<<"ANASEN Geometry File -- Coordinates for Detectors"<<std::endl;
	output<<"Edges: x y z"<<std::endl;
	for(unsigned int i=0; i<x.size(); i++) {
		output<<x[i]<<" "<<y[i]<<" "<<z[i]<<std::endl;
	}
	output<<"Centers: x y z"<<std::endl;
	for(unsigned int i=0; i<cx.size(); i++) {
		output<<cx[i]<<" "<<cy[i]<<" "<<cz[i]<<std::endl;
	}

	output.close();
}

double AnasenEfficiency::RunConsistencyCheck() {
	std::vector<Mask::Vec3> r1_points;
	std::vector<Mask::Vec3> r2_points;
	std::vector<Mask::Vec3> fqqq_points;
	std::vector<Mask::Vec3> bqqq_points;
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			r1_points.push_back(m_Ring1[i].GetHitCoordinates(j, 0));
		}
	}
	for(int i=0; i<n_sx3_per_ring; i++) {
		for(int j=0; j<4; j++) {
			r2_points.push_back(m_Ring2[i].GetHitCoordinates(j, 0));
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<16; k++) {
				fqqq_points.push_back(m_forwardQQQs[i].GetHitCoordinates(j, k));
			}
		}
	}
	for(int i=0; i<n_qqq; i++) {
		for(int j=0; j<16; j++) {
			for(int k=0; k<16; k++) {
				bqqq_points.push_back(m_backwardQQQs[i].GetHitCoordinates(j, k));
			}
		}
	}

	int npoints = r1_points.size() + r2_points.size() + fqqq_points.size() + bqqq_points.size();
	int count = 0;
	Mask::Vec3 coords;
	for(auto& point : r1_points) {
		for(auto& sx3 : m_Ring1) {
			auto result = sx3.GetChannelRatio(point.GetTheta(), point.GetPhi());
			coords = sx3.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}
	for(auto& point : r2_points) {
		for(auto& sx3 : m_Ring2) {
			auto result = sx3.GetChannelRatio(point.GetTheta(), point.GetPhi());
			coords = sx3.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}
	for(auto& point : fqqq_points) {
		for(auto& qqq : m_forwardQQQs) {
			auto result = qqq.GetTrajectoryRingWedge(point.GetTheta(), point.GetPhi());
			coords = qqq.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}
	for(auto& point : bqqq_points) {
		for(auto& qqq : m_backwardQQQs) {
			auto result = qqq.GetTrajectoryRingWedge(point.GetTheta(), point.GetPhi());
			coords = qqq.GetHitCoordinates(result.first, result.second);
			if(IsDoubleEqual(point.GetX(), coords.GetX()) && IsDoubleEqual(point.GetY(), coords.GetY()) && IsDoubleEqual(point.GetZ(), coords.GetZ())) {
				count++;
				break;
			}
		}
	}

	double ratio = ((double)count)/((double)npoints);

	return ratio;

}

bool AnasenEfficiency::IsRing1(double theta, double phi) {
	for(auto& sx3 : m_Ring1) {
		auto result = sx3.GetChannelRatio(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}

	return false;
}

bool AnasenEfficiency::IsRing2(double theta, double phi) {
	for(auto& sx3 : m_Ring2) {
		auto result = sx3.GetChannelRatio(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}

	return false;
}

bool AnasenEfficiency::IsQQQ(double theta, double phi) {
	for(auto& qqq : m_forwardQQQs) {
		auto result = qqq.GetTrajectoryRingWedge(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}

	
	for(auto& qqq : m_backwardQQQs) {
		auto result = qqq.GetTrajectoryRingWedge(theta, phi);
		if(result.first != -1) {
			return true;
		}
	}



	return false;
}

void AnasenEfficiency::CalculateEfficiency(const std::string& inputname, const std::string& outputname, const std::string& statsname) {
	std::cout<<"----------ANASEN Efficiency Calculation----------"<<std::endl;
	std::cout<<"Loading in output from kinematics simulation: "<<inputname<<std::endl;
	std::cout<<"Running efficiency calculation..."<<std::endl;

	Mask::MaskFile input(inputname, Mask::MaskFile::FileType::read);
	Mask::MaskFile output(outputname, Mask::MaskFile::FileType::write);
	std::ofstream stats(statsname);
	stats<<std::setprecision(5);

	Mask::MaskFileHeader header = input.ReadHeader();
	output.WriteHeader(header.rxn_type, header.nsamples);
	stats<<"Efficiency statistics for data from "<<inputname<<" using the ANASEN geometry"<<std::endl;
	stats<<"Given in order of target=0, projectile=1, ejectile=2, residual=3, .... etc."<<std::endl;

	Mask::MaskFileData data;

	std::vector<int> counts;
	switch(header.rxn_type) {
		case 0:
			counts.resize(3, 0);
			break;
		case 1:
			counts.resize(4, 0);
			break;
		case 2:
			counts.resize(6, 0);
			break;
		case 3:
			counts.resize(8, 0);
			break;
		default:
		{
			std::cerr<<"Bad reaction type at AnasenEfficiency::CalculateEfficiency (given value: "<<header.rxn_type<<"). Quiting..."<<std::endl;
			input.Close();
			output.Close();
			stats.close();
			return;
		}
	}

	int percent5 = header.nsamples*0.05;
	int count = 0;
	int npercent = 0;

	while(true) {
		if(++count == percent5) {//Show progress every 5%
			npercent++;
			count = 0;
			std::cout<<"\rPercent completed: "<<npercent*5<<"%"<<std::flush;
		}

		data = input.ReadData();
		if(data.eof)
			break;

		for(unsigned int i=0; i<data.Z.size(); i++) {
			if(data.KE[i] >= threshold && (IsRing1(data.theta[i], data.phi[i]) || IsRing2(data.theta[i], data.phi[i]) || IsQQQ(data.theta[i], data.phi[i]))) {
				data.detect_flag[i] = true;
				counts[i]++;
			} else if(data.detect_flag[i] == true) {
				data.detect_flag[i] = false;
			}
		}

		output.WriteData(data);
	}

	input.Close();
	output.Close();

	stats<<std::setw(10)<<"Index"<<"|"<<std::setw(10)<<"Efficiency"<<std::endl;
	stats<<"---------------------"<<std::endl;
	for(unsigned int i=0; i<counts.size(); i++) {
		stats<<std::setw(10)<<i<<"|"<<std::setw(10)<<((double)counts[i]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	stats.close();

	std::cout<<std::endl;
	std::cout<<"Complete."<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
}