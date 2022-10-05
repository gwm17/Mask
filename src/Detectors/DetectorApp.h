#ifndef DETECTOR_APP_H
#define DETECTOR_APP_H

#include "DetectorArray.h"
#include "Mask/FileWriter.h"
#include "Mask/FileReader.h"
#include "Mask/ThreadPool.h"

#include <string>
#include <vector>
#include <memory>

class DetectorApp
{
public:
    DetectorApp();
    ~DetectorApp();

    bool LoadConfig(const std::string& filename);
    
    void Run();

private:

    std::vector<DetectorArray*> m_detectorList; //One array per thread
    std::vector<uint64_t> m_chunkSamples;
    Mask::FileWriter m_fileWriter;
    Mask::FileReader m_fileReader;

    std::string m_inputFileName;
    std::string m_outputFileName;
    std::string m_deadChannelFileName;

    uint64_t m_nthreads;
    uint64_t m_nentries;

    std::unique_ptr<Mask::ThreadPool<DetectorArray*, uint64_t>> m_resources;
};

#endif