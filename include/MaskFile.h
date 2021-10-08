#ifndef MASKFILE_H
#define MASKFILE_H

#include <string>
#include <fstream>
#include <vector>

#include "Nucleus.h"
#include "RxnType.h"

namespace Mask {

	struct MaskFileHeader {
		RxnType rxn_type = RxnType::None;
		int nsamples = -1;
	};

	struct MaskFileData {
		std::vector<double> E, KE, p, theta, phi; //ordered: target, (if not decay)projectile, ejectile, residual, break1...
		std::vector<int> Z, A;
		std::vector<bool> detect_flag;
		bool eof = false; //flag on end of file
	};

	class MaskFile {
	public:
		enum class FileType {
			read,
			write,
			append,
			none
		};
	
		MaskFile();
		MaskFile(const std::string& name, MaskFile::FileType type);
		bool Open(const std::string& name, MaskFile::FileType type);
		inline bool IsOpen() { return file.is_open(); }
		void Close();
		
		void WriteHeader(RxnType rxn_type, int nsamples);
		void WriteData(std::vector<Nucleus>& data);
		void WriteData(MaskFileData& data);
		MaskFileHeader ReadHeader();
		MaskFileData ReadData();
	
		
	
	private:

		FileType file_type;
		std::string filename;
		uint32_t buffer_position;
		uint32_t buffer_end;
		uint32_t data_size;
		RxnType m_rxn_type;
		uint32_t buffersize_bytes;
	
		std::fstream file;
	
		std::vector<char> data_buffer;
	
		static constexpr uint32_t onestep_rxn_n = 2;
		static constexpr uint32_t twostep_rxn_n = 4;
		static constexpr uint32_t threestep_rxn_n = 6;
	
		static constexpr uint64_t buffersize = 10000; //number of events
		static constexpr int width = 0;
		static constexpr int precision = 3;

		static constexpr std::size_t double_size = sizeof(double);
		static constexpr std::size_t int_size = sizeof(uint32_t);
		static constexpr std::size_t bool_size = sizeof(bool);
	};

};

#endif