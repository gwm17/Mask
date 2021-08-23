#include "MaskFile.h"
#include <iostream>
#include <iomanip>
#include <sstream>

namespace Mask {

MaskFile::MaskFile() : 
	file_type(MaskFile::FileType::none), filename(""), buffer_position(0), count(0), file()
{
}

MaskFile::MaskFile(const std::string& name, MaskFile::FileType type) :
	file_type(type), filename(name), buffer_position(0), count(0), file()
{
	Open(filename, type);
}

bool MaskFile::Open(const std::string& name, MaskFile::FileType type) {
	if(IsOpen()) {
		std::cerr<<"Attempted to open file that is already open!"<<std::endl;
		return false;
	}

	if(type == FileType::read) {
		file.open(name, std::ios::in);
		buffer_position = 0;
		if(IsOpen()) ReadHeader();
	} else if (type == FileType::write) {
		file.open(name, std::ios::out | std::ios::trunc);
		file<<std::setprecision(precision);
		file<<std::scientific;
	} else if (type == FileType::append) {
		file.open(name, std::ios::out | std::ios::app);
		file<<std::setprecision(precision);
		file<<std::scientific;
	} else {
		std::cerr<<"Invalid FileType at MaskFile::Open()"<<std::endl;
		return IsOpen();
	}

	count = 0;

	return IsOpen();
}

void MaskFile::Close() {
	if(data_buffer.size() > 0) {
		FlushBuffer();
	}
	file.close();
}

void MaskFile::WriteHeader(std::vector<Nucleus>& nuclei) {
	file<<std::setw(width)<<"ZT"<<","
		<<std::setw(width)<<"AT"<<","
		<<std::setw(width)<<"ZP"<<","
		<<std::setw(width)<<"AP"<<","
		<<std::setw(width)<<"ZE"<<","
		<<std::setw(width)<<"AE"<<","
		<<std::setw(width)<<"ZR"<<","
		<<std::setw(width)<<"AR"<<","
		<<std::setw(width)<<"ZB1"<<","
		<<std::setw(width)<<"AB1"<<","
		<<std::setw(width)<<"ZB2"<<","
		<<std::setw(width)<<"AB2"<<","
		<<std::setw(width)<<"ZB3"<<","
		<<std::setw(width)<<"AB3"<<","
		<<std::setw(width)<<"ZB4"<<","
		<<std::setw(width)<<"AB4"<<","
		<<std::setw(width)<<"KEP"<<","
		<<std::endl;

	file<<std::setw(width)<<nuclei[0].GetZ()<<","
		<<std::setw(width)<<nuclei[0].GetA()<<","
		<<std::setw(width)<<nuclei[1].GetZ()<<","
		<<std::setw(width)<<nuclei[1].GetA()<<","
		<<std::setw(width)<<nuclei[2].GetZ()<<","
		<<std::setw(width)<<nuclei[2].GetA()<<","
		<<std::setw(width)<<nuclei[3].GetZ()<<","
		<<std::setw(width)<<nuclei[3].GetA()<<",";
	if(nuclei.size() == 6) {
		file<<std::setw(width)<<nuclei[4].GetZ()<<","
			<<std::setw(width)<<nuclei[4].GetA()<<","
			<<std::setw(width)<<nuclei[5].GetZ()<<","
			<<std::setw(width)<<nuclei[5].GetA()<<",";
	} else if(nuclei.size() == 8) {
		file<<std::setw(width)<<nuclei[4].GetZ()<<","
			<<std::setw(width)<<nuclei[4].GetA()<<","
			<<std::setw(width)<<nuclei[5].GetZ()<<","
			<<std::setw(width)<<nuclei[5].GetA()<<","
			<<std::setw(width)<<nuclei[6].GetZ()<<","
			<<std::setw(width)<<nuclei[6].GetA()<<","
			<<std::setw(width)<<nuclei[7].GetZ()<<","
			<<std::setw(width)<<nuclei[7].GetA()<<",";
	}
	file<<std::setw(width)<<nuclei[1].GetKE()<<",";
	file<<std::endl;

	file<<std::setw(width)<<"Event"<<","
		<<std::setw(width)<<"EE"<<","<<std::setw(width)<<"KEE"<<","<<std::setw(width)<<"PE"<<","<<std::setw(width)<<"ThetaE"<<","<<std::setw(width)<<"PhiE"<<","<<std::setw(width)<<"ThetaCME"<<","
		<<std::setw(width)<<"ER"<<","<<std::setw(width)<<"KER"<<","<<std::setw(width)<<"PR"<<","<<std::setw(width)<<"ThetaR"<<","<<std::setw(width)<<"PhiR"<<","<<std::setw(width)<<"ThetaCMR"<<","
		<<std::setw(width)<<"EB1"<<","<<std::setw(width)<<"KEB1"<<","<<std::setw(width)<<"PB1"<<","<<std::setw(width)<<"ThetaB1"<<","<<std::setw(width)<<"PhiB1"<<","<<std::setw(width)<<"ThetaCMB1"<<","
		<<std::setw(width)<<"EB2"<<","<<std::setw(width)<<"KEB2"<<","<<std::setw(width)<<"PB2"<<","<<std::setw(width)<<"ThetaB2"<<","<<std::setw(width)<<"PhiB2"<<","<<std::setw(width)<<"ThetaCMB2"<<","
		<<std::setw(width)<<"EB3"<<","<<std::setw(width)<<"KEB3"<<","<<std::setw(width)<<"PB3"<<","<<std::setw(width)<<"ThetaB3"<<","<<std::setw(width)<<"PhiB3"<<","<<std::setw(width)<<"ThetaCMB3"<<","
		<<std::setw(width)<<"EB4"<<","<<std::setw(width)<<"KEB4"<<","<<std::setw(width)<<"PB4"<<","<<std::setw(width)<<"ThetaB4"<<","<<std::setw(width)<<"PhiB4"<<","<<std::setw(width)<<"ThetaCMB4"<<","
		<<std::endl;
}

void MaskFile::ReadHeader() {
	std::string junk;
	std::getline(file, junk);
}

void MaskFile::WriteData(std::vector<Nucleus>& data) {
	if(count == 0 && data_buffer.size() == 0) {
		WriteHeader(data);
	}
	
	data_buffer.push_back(data);
	if(data_buffer.size() == buffersize) {
		FlushBuffer();
	}
}

void MaskFile::FlushBuffer() {
	for(auto& event : data_buffer) {
		file<<std::setw(width)<<count<<",";
		for(unsigned int i=2; i< event.size(); i++) {
			Nucleus& nuc = event[i];
			file<<std::setw(width)<<nuc.GetE()<<","
				<<std::setw(width)<<nuc.GetKE()<<","
				<<std::setw(width)<<nuc.GetP()<<","
				<<std::setw(width)<<nuc.GetTheta()<<","
				<<std::setw(width)<<nuc.GetPhi()<<","
				<<std::setw(width)<<nuc.GetThetaCM()<<",";
		}
		file<<std::endl;
		count++;
	}
	data_buffer.clear();
}

std::vector<Nucleus> MaskFile::ReadData() {
	if(buffer_position == data_buffer.size()) {
		FillBuffer();
		buffer_position = 0;
	}
	std::vector<Nucleus> data = data_buffer[buffer_position];
	buffer_position++;
	
	return data;
}

void MaskFile::FillBuffer() {
	std::vector<Nucleus> data;
	std::string line;
	std::stringstream linebuf;
	std::string junk;
	int Z, A;
	double E, p, theta, phi, tcm;
	data_buffer.clear();

	while(data_buffer.size() <= buffersize) {
		if(!std::getline(file, line)) break;
		linebuf.str(line);
		while(linebuf>>Z) {
			linebuf>>A;
			linebuf>>E>>junk>>junk>>p>>theta>>phi>>tcm;
			data.emplace_back(Z, A);
			data[data.size()-1].SetVectorSpherical(theta, phi, p, E);
			data[data.size()-1].SetThetaCM(tcm);
		}
		data_buffer.push_back(data);
		data.clear();
	}
}

};