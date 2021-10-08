#include "AnasenEfficiency.h"
#include "Kinematics.h"
#include <fstream>
#include <iomanip>

AnasenEfficiency::AnasenEfficiency() :
	DetectorEfficiency(), det_silicon(si_thickness)
{
	for(int i=0; i<n_sx3_per_ring; i++) {
		m_Ring1.emplace_back(4, sx3_length, sx3_width, ring_phi[i], ring1_z, ring_rho[i]);
		m_Ring1[i].TurnOnRandomizedCoordinates();
		m_Ring2.emplace_back(4, sx3_length, sx3_width, ring_phi[i], ring2_z, ring_rho[i]);
		m_Ring2[i].TurnOnRandomizedCoordinates();
	}
	for(int i=0; i<n_qqq; i++) {
		m_forwardQQQs.emplace_back(qqq_rinner, qqq_router, qqq_deltaphi, qqq_phi[i], qqq_z[i]);
		m_forwardQQQs[i].TurnOnRandomizedCoordinates();
		m_backwardQQQs.emplace_back(qqq_rinner, qqq_router, qqq_deltaphi, qqq_phi[i], (-1.0)*qqq_z[i]);
		m_backwardQQQs[i].TurnOnRandomizedCoordinates();
	}

	std::vector<int> det_z = {14};
	std::vector<int> det_a = {28};
	std::vector<int> det_stoich = {1};
	det_silicon.SetElements(det_z, det_a, det_stoich);

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
	for(unsigned int i=0; i<x.size(); i++) {
		output<<x[i]<<" "<<y[i]<<" "<<z[i]<<std::endl;
	}
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

DetectorResult AnasenEfficiency::IsRing1(Mask::Nucleus& nucleus) {

	DetectorResult observation;
	//Mask::Vec3 coords;
	double thetaIncident, eloss, e_dep;
	for(auto& sx3 : m_Ring1) {
		auto result = sx3.GetChannelRatio(nucleus.GetTheta(), nucleus.GetPhi());
		if(result.first != -1) {
			//coords = sx3.GetHitCoordinates(result.first, result.second);
			observation.detectFlag = true;
			observation.direction = sx3.GetHitCoordinates(result.first, result.second);
			thetaIncident = std::acos(observation.direction.Dot(sx3.GetNormRotated())/observation.direction.GetR());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			//e_dep = det_silicon.getEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			observation.energy_deposited = det_silicon.GetEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			observation.det_name = "R1";
			//return std::make_pair(true, e_dep);
			return observation;
		}
	}

	//return std::make_pair(false, 0.0);
	return observation;
}

DetectorResult AnasenEfficiency::IsRing2(Mask::Nucleus& nucleus) {

	DetectorResult observation;
	//Mask::Vec3 coords;
	double thetaIncident, eloss, e_dep;
	for(auto& sx3 : m_Ring2) {
		auto result = sx3.GetChannelRatio(nucleus.GetTheta(), nucleus.GetPhi());
		if(result.first != -1) {
			//coords = sx3.GetHitCoordinates(result.first, result.second);
			observation.detectFlag = true;
			observation.direction = sx3.GetHitCoordinates(result.first, result.second);
			thetaIncident = std::acos(observation.direction.Dot(sx3.GetNormRotated())/observation.direction.GetR());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			//e_dep = det_silicon.getEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			observation.energy_deposited = det_silicon.GetEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			observation.det_name = "R2";
			//return std::make_pair(true, e_dep);
			return observation;
		}
	}

	return observation;
}

DetectorResult AnasenEfficiency::IsQQQ(Mask::Nucleus& nucleus) {

	DetectorResult observation;
	//Mask::Vec3 coords;
	double thetaIncident, eloss, e_dep;

	for(auto& qqq : m_forwardQQQs) {
		auto result = qqq.GetTrajectoryRingWedge(nucleus.GetTheta(), nucleus.GetPhi());
		if(result.first != -1) {
			//coords = qqq.GetHitCoordinates(result.first, result.second);
			observation.detectFlag = true;
			observation.direction = qqq.GetHitCoordinates(result.first, result.second);
			thetaIncident = std::acos(observation.direction.Dot(qqq.GetNorm())/observation.direction.GetR());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			//e_dep = det_silicon.getEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			observation.energy_deposited = det_silicon.GetEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			//return std::make_pair(true, e_dep);
			observation.det_name = "FQQQ";
			return observation;
		}
	}

	
	for(auto& qqq : m_backwardQQQs) {
		auto result = qqq.GetTrajectoryRingWedge(nucleus.GetTheta(), nucleus.GetPhi());
		if(result.first != -1) {
			//coords = qqq.GetHitCoordinates(result.first, result.second);
			observation.detectFlag = true;
			observation.direction = qqq.GetHitCoordinates(result.first, result.second);
			thetaIncident = std::acos(observation.direction.Dot(qqq.GetNorm())/observation.direction.GetR());
			if(thetaIncident > M_PI/2.0)
				thetaIncident = M_PI - thetaIncident;

			//e_dep = det_silicon.getEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			observation.energy_deposited = det_silicon.GetEnergyLossTotal(nucleus.GetZ(), nucleus.GetA(), nucleus.GetKE(), thetaIncident);
			//return std::make_pair(true, e_dep);
			observation.det_name = "BQQQ";
			return observation;
		}
	}

	//return std::make_pair(false, 0.0);
	return observation;
}

DetectorResult AnasenEfficiency::IsAnasen(Mask::Nucleus& nucleus) {
	DetectorResult result;
	if(nucleus.GetKE() <= threshold)
		return result;

	if(!result.detectFlag) {
		result = IsRing1(nucleus);
	}
	if(!result.detectFlag) {
		result = IsRing2(nucleus);
	}
	if(!result.detectFlag) {
		result = IsQQQ(nucleus);
	}

	return result;
}

void AnasenEfficiency::CountCoincidences(const Mask::MaskFileData& data, std::vector<int>& counts, int rxn_type) {
	if (rxn_type == 0 && data.detect_flag[1] && data.detect_flag[2])
	{
		counts[0]++;
	}
	else if (rxn_type == 1 && data.detect_flag[2] && data.detect_flag[3])
	{
		counts[0]++;
	}
	else if(rxn_type == 2)
	{
		if(data.detect_flag[2] && data.detect_flag[4]) 
		{
			counts[0]++;
		}
		if(data.detect_flag[2] && data.detect_flag[5])
		{
			counts[1]++;
		}
		if(data.detect_flag[4] && data.detect_flag[5])
		{
			counts[2]++;
		}
		if(data.detect_flag[2] && data.detect_flag[4] && data.detect_flag[5])
		{
			counts[3]++;
		}
	}
	else if(rxn_type == 3)
	{
		if(data.detect_flag[2] && data.detect_flag[4]) 
		{
			counts[0]++;
		}
		if(data.detect_flag[2] && data.detect_flag[6])
		{
			counts[1]++;
		}
		if(data.detect_flag[2] && data.detect_flag[7])
		{
			counts[2]++;
		}
		if(data.detect_flag[4] && data.detect_flag[6])
		{
			counts[3]++;
		}
		if(data.detect_flag[4] && data.detect_flag[7])
		{
			counts[4]++;
		}
		if(data.detect_flag[6] && data.detect_flag[7])
		{
			counts[5]++;
		}
		if(data.detect_flag[2] && data.detect_flag[4] && data.detect_flag[6])
		{
			counts[6]++;
		}
		if(data.detect_flag[2] && data.detect_flag[4] && data.detect_flag[7])
		{
			counts[7]++;
		}
		if(data.detect_flag[2] && data.detect_flag[6] && data.detect_flag[7])
		{
			counts[8]++;
		}
		if(data.detect_flag[4] && data.detect_flag[6] && data.detect_flag[7])
		{
			counts[9]++;
		}
		if(data.detect_flag[2] && data.detect_flag[4] && data.detect_flag[6] && data.detect_flag[7])
		{
			counts[10]++;
		}
	}
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
	std::vector<int> coinc_counts;
	switch(header.rxn_type) {
		case 0:
			counts.resize(3, 0);
			coinc_counts.resize(1, 0);
			break;
		case 1:
			counts.resize(4, 0);
			coinc_counts.resize(1, 0);
			break;
		case 2:
			counts.resize(6, 0);
			coinc_counts.resize(4, 0);
			break;
		case 3:
			counts.resize(8, 0);
			coinc_counts.resize(11, 0);
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

	Mask::Nucleus nucleus;
	int index=0;
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
			nucleus.SetIsotope(data.Z[i], data.A[i]);
			nucleus.SetVectorSpherical(data.theta[i], data.phi[i], data.p[i], data.E[i]);
			auto result = IsAnasen(nucleus);
			if(result.detectFlag) {
				data.detect_flag[i] = true;
				data.KE[i] = result.energy_deposited;
				data.theta[i] = result.direction.GetTheta();
				data.phi[i] = result.direction.GetPhi();
				counts[i]++;
			} else if(data.detect_flag[i] == true) {
				data.detect_flag[i] = false;
			}
		}

		CountCoincidences(data, coinc_counts, header.rxn_type);

		output.WriteData(data);

		index++;
	}

	input.Close();
	output.Close();

	stats<<std::setw(10)<<"Index"<<"|"<<std::setw(10)<<"Efficiency"<<std::endl;
	stats<<"---------------------"<<std::endl;
	for(unsigned int i=0; i<counts.size(); i++) {
		stats<<std::setw(10)<<i<<"|"<<std::setw(10)<<((double)counts[i]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	stats<<"Coincidence Efficiency"<<std::endl;
	stats<<"---------------------"<<std::endl;
	if(header.rxn_type == 0)
	{
		stats<<std::setw(10)<<"1 + 2"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(header.rxn_type == 1)
	{
		stats<<std::setw(10)<<"2 + 3"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(header.rxn_type == 2)
	{
		stats<<std::setw(10)<<"2 + 4"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 5"<<"|"<<std::setw(10)<<((double)coinc_counts[1]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 5"<<"|"<<std::setw(10)<<((double)coinc_counts[2]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 5"<<"|"<<std::setw(10)<<((double)coinc_counts[3]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	else if(header.rxn_type == 3)
	{
		stats<<std::setw(10)<<"2 + 4"<<"|"<<std::setw(10)<<((double)coinc_counts[0]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 6"<<"|"<<std::setw(10)<<((double)coinc_counts[1]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[2]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 6"<<"|"<<std::setw(10)<<((double)coinc_counts[3]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[4]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[5]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 6"<<"|"<<std::setw(10)<<((double)coinc_counts[6]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[7]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[8]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"4 + 6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[9]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
		stats<<std::setw(10)<<"2 + 4 + 6 + 7"<<"|"<<std::setw(10)<<((double)coinc_counts[10]/header.nsamples)<<std::endl;
		stats<<"---------------------"<<std::endl;
	}
	stats.close();

	std::cout<<std::endl;
	std::cout<<"Complete."<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
}