#ifndef MASKFILE_H
#define MASKFILE_H

#include <string>
#include <fstream>
#include <vector>

#include "Nucleus.h"

namespace Mask {

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
	
	void WriteHeader(std::vector<Nucleus>& data);
	void WriteData(std::vector<Nucleus>& data);
	void ReadHeader();
	std::vector<Nucleus> ReadData();

	

private:
	void FlushBuffer();
	void FillBuffer();

	FileType file_type;
	std::string filename;
	unsigned int buffer_position;
	int count;

	std::fstream file;

	std::vector<std::vector<Nucleus>> data_buffer;

	static constexpr int buffersize = 1000;
	static constexpr int width = 0;
	static constexpr int precision = 3;
};

};

#endif