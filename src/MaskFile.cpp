#include "MaskFile.h"
#include <iostream>
#include <iomanip>
#include <sstream>

/*

	FORMAT

	HEADER (contains rxntype & nuclei numbers & beam kinetic energy) (64 bits, 8 bytes)
		NSAMPLES(32bit) RXNTYPE(32bit)
	end HEADER

	There are NSAMPLES * (number of saved nuclei) data in the file. The number of nuclei saved is related to the 
	RXNTYPE. All nuclei (target, projectile, ejectile, residual, break1, etc...) are saved. A datum is as follows:

	DATUM (contains kinematic data for a nucleus) (384 bits, 48 bytes)
		Z(32bit) A(32bit) DETECTFLAG(32bit) E(64bit) KE(64bit) P(64bit) THETA(64bit) PHI(64bit)
	end DATUM


*/

namespace Mask {

	MaskFile::MaskFile() : 
		file_type(MaskFile::FileType::none), filename(""), buffer_position(0), buffer_end(0), data_size(0), buffersize_bytes(0), file()
	{
	}
	
	MaskFile::MaskFile(const std::string& name, MaskFile::FileType type) :
		file_type(type), filename(name), buffer_position(0), buffer_end(0), data_size(0), buffersize_bytes(0), file()
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
		} else if (type == FileType::write) {
			file.open(name, std::ios::out | std::ios::trunc);
		} else if (type == FileType::append) {
			file.open(name, std::ios::out | std::ios::app);
		} else {
			std::cerr<<"Invalid FileType at MaskFile::Open()"<<std::endl;
			return IsOpen();
		}
	
	
		return IsOpen();
	}
	
	void MaskFile::Close() {
		//Final flush (if necessary)
		if(buffer_position > 0 && buffer_position < buffersize_bytes && (file_type == FileType::write || file_type == FileType::append)) {
			file.write(data_buffer.data(), buffer_position);
		}
		file.close();
	}
	
	void MaskFile::WriteHeader(int rxn_type, int nsamples) {
		//size of a datum = data nuclei per event * (# doubles per nucleus * sizeof double + # ints per nucleus * sizeof int + sizeof bool)
		if(rxn_type == 0) {
			m_rxn_type = 0;
			data_size = 3 * ( 5 * double_size + 2 * int_size + bool_size); 
		}
		else if (rxn_type == 1) {
			m_rxn_type = 1;
			data_size = 4 * ( 5 * double_size + 2 * int_size + bool_size);
		}
		else if (rxn_type == 2) {
			m_rxn_type = 2;
			data_size = 6 * ( 5 * double_size + 2 * int_size + bool_size);
		}
		else if (rxn_type == 3) {
			m_rxn_type = 3;
			data_size = 8 * ( 5 * double_size + 2 * int_size + bool_size);
		} else {
			std::cerr<<"Invalid number of nuclei at MaskFile::WriteHeader! Returning"<<std::endl;
			return;
		}
		buffersize_bytes = buffersize * data_size; //buffersize_bytes = # of events * size of an event
	
		data_buffer.resize(buffersize_bytes);
	
		file.write((char*) &nsamples, int_size);
		file.write((char*) &m_rxn_type, int_size);
	}
	
	MaskFileHeader MaskFile::ReadHeader() {
		MaskFileHeader header;
		std::vector<char> temp_buffer(4);
		file.read(temp_buffer.data(), 4);
		header.nsamples = *(int*)(&temp_buffer[0]);
		file.read(temp_buffer.data(), 4);
		m_rxn_type = *(int*)(&temp_buffer[0]);
	
		//size of a datum = data nuclei per event * (# doubles per nucleus * sizeof double + # ints per nucleus * sizeof int + sizeof bool)
		if(m_rxn_type == 0) {
			data_size = 3 * ( 5 * double_size + 2 * int_size + bool_size);
		}
		else if (m_rxn_type == 1) {
			data_size = 4 * ( 5 * double_size + 2 * int_size + bool_size);
		}
		else if (m_rxn_type == 2) {
			data_size = 6 * ( 5 * double_size + 2 * int_size + bool_size);
		}
		else if (m_rxn_type == 3) {
			data_size = 8 * ( 5 * double_size + 2 * int_size + bool_size);
		} else {
			std::cerr<<"Unexpected reaction type at MaskFile::ReadHeader (rxn type = "<<m_rxn_type<<")! Returning"<<std::endl;
			return header;
		}
		buffersize_bytes = buffersize * data_size;//buffersize_bytes = size of a datum * # of events
	
		header.rxn_type = m_rxn_type;
	
		data_buffer.resize(buffersize_bytes);
	
		return header;
	}
	
	void MaskFile::WriteData(std::vector<Nucleus>& data) {
	
		char* data_pointer;
		double datum;
		int number;
		bool flag;
		std::size_t j;
		for(unsigned int i=0; i<data.size(); i++) {
			number = data[i].GetZ();
			data_pointer = (char*) &number;
			for(j=0; j<int_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			number = data[i].GetA();
			data_pointer = (char*) &number;
			for(j=0; j<int_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			flag = data[i].IsDetected();
			data_pointer = (char*) &flag;
			for(j=0; j<bool_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			datum = data[i].GetE();
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			datum = data[i].GetKE();
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
	
			datum = data[i].GetP();
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
	
			datum = data[i].GetTheta();
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
	
			datum = data[i].GetPhi();
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
		}
	
		//Flush the buffer when it is full, and reset the position.
		if(buffer_position == buffersize_bytes) {
			file.write(data_buffer.data(), data_buffer.size());
			buffer_position = 0;
		}
	}

	void MaskFile::WriteData(MaskFileData& data) {
	
		char* data_pointer;
		double datum;
		int number;
		bool flag;
		std::size_t j;
		for(unsigned int i=0; i<data.Z.size(); i++) {
			number = data.Z[i];
			data_pointer = (char*) &number;
			for(j=0; j<int_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			number = data.A[i];
			data_pointer = (char*) &number;
			for(j=0; j<int_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			flag = data.detect_flag[i];
			data_pointer = (char*) &flag;
			for(j=0; j<bool_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			datum = data.E[i];
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}

			datum = data.KE[i];
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
	
			datum = data.p[i];
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
	
			datum = data.theta[i];
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
	
			datum = data.phi[i];
			data_pointer = (char*) &datum;
			for(j=0; j<double_size; j++) {
				data_buffer[buffer_position] = *(data_pointer + j);
				buffer_position++;
			}
		}
	
		//Flush the buffer when it is full, and reset the position.
		if(buffer_position == buffersize_bytes) {
			file.write(data_buffer.data(), data_buffer.size());
			buffer_position = 0;
		}
	}
	
	/*
		Read data from the buffer and submit it to the client side as a MaskFileData struct.
		When file reaches the end of the file (no more data to read), an empty MaskFileData with 
		eof == true is sent out signaling that the file is finished.

		Should be used like

		Mask::MaskFile input(file, Mask::MaskFile::FileType::read);
		Mask::MaskFileHeader header = input.ReadHeader();
		Mask::MaskFileData data;
		while(true) {
			data = input.ReadData();
			if(data.eof) break;

			Do some stuff...

		}
		input.Close();

		Good luck
	*/
	MaskFileData MaskFile::ReadData() {
	
		MaskFileData data;
	
		//Fill the buffer when needed, reset the positon, and find the end
		if(buffer_position == data_buffer.size() || buffer_position == buffer_end) {
			file.read(data_buffer.data(), buffersize_bytes);
			buffer_position = 0;
			buffer_end = file.gcount();
			if(buffer_end == 0 && file.eof()) {
				data.eof = true;
				return data;
			}
		}
	
		unsigned int local_end = buffer_position + data_size;
		if(local_end > buffer_end) {
			std::cerr<<"Attempting to read past end of file at MaskFile::ReadData! Returning empty"<<std::endl;
			data.eof = true;
			return data;
		}

		while(buffer_position < local_end) {
			data.Z.push_back(*(int*)(&data_buffer[buffer_position]));
			buffer_position += int_size;
			data.A.push_back(*(int*)(&data_buffer[buffer_position]));
			buffer_position += int_size;
			data.detect_flag.push_back(*(bool*)(&data_buffer[buffer_position]));
			buffer_position += bool_size;
			data.E.push_back(*(double*)(&data_buffer[buffer_position]));
			buffer_position += double_size;
			data.KE.push_back(*(double*)(&data_buffer[buffer_position]));
			buffer_position += double_size;
			data.p.push_back(*(double*)(&data_buffer[buffer_position]));
			buffer_position += double_size;
			data.theta.push_back(*(double*)(&data_buffer[buffer_position]));
			buffer_position += double_size;
			data.phi.push_back(*(double*)(&data_buffer[buffer_position]));
			buffer_position += double_size;
		}
	
		return data;
	}

};