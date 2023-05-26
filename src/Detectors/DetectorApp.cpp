#include "DetectorApp.h"
#include <fstream>
#include <iostream>

#include "yaml-cpp/yaml.h"

DetectorApp::DetectorApp() :
    m_resources(nullptr)
{
}

DetectorApp::~DetectorApp()
{
    for(std::size_t i=0; i<m_detectorList.size(); i++)
        delete m_detectorList[i];
}

bool DetectorApp::ParseConfig(const std::string& filename)
{
    std::cout<<"----------Detector Efficiency Calculation----------"<<std::endl;
    YAML::Node data;
    try
    {
        data = YAML::LoadFile(filename);
    }
    catch (YAML::ParserException& e)
    {
        std::cerr << "Could not load config file " << filename  << " with error: " << e.what() << std::endl;
        return false;
    }

    m_inputFileName = data["InputDataFile"].as<std::string>();
    m_outputFileName = data["OutputDataFile"].as<std::string>();
    m_deadChannelFileName = data["DeadChannelFile"].as<std::string>();
    m_nthreads = data["NumberOfThreads"].as<uint64_t>();
    ArrayType type = StringToArrayType(data["ArrayType"].as<std::string>());

    std::cout << "Creating " << m_nthreads << " detector arrays..." << std::endl;
    for(uint64_t i=0; i<m_nthreads; i++)
    {
        m_detectorList.push_back(CreateDetectorArray(type));
        if(m_deadChannelFileName != "None")
            m_detectorList.back()->SetDeadChannelMap(m_deadChannelFileName);
    }

    std::cout << "Allocating " << m_nthreads << " threads..." << std::endl;
    m_resources = std::make_unique<Mask::ThreadPool<DetectorArray*, uint64_t>>(m_nthreads);

    std::cout << "Opening input data file " << m_inputFileName << "..." << std::endl;
    m_fileReader.Open(m_inputFileName, "SimTree");
    if(!m_fileReader.IsOpen() || !m_fileReader.IsTree())
    {
        std::cerr << "Unable to open input data file " << m_inputFileName << std::endl;
        return false;
    }
    m_nentries = m_fileReader.GetSize();
    std::cout << "Detected " << m_nentries << " events in the file..." << std::endl;

    //Little bit of integer division mangling to make sure we read every event in file
    uint64_t quotient = m_nentries / m_nthreads;
    uint64_t remainder = m_nentries % m_nthreads;
    m_chunkSamples.push_back(quotient + remainder);
    for(uint64_t i=1; i<m_nthreads; i++)
        m_chunkSamples.push_back(quotient);

    std::cout << "Opening output data file " << m_outputFileName << std::endl;
    m_fileWriter.Open(m_outputFileName, "SimTree");
    if(!m_fileWriter.IsOpen() || !m_fileWriter.IsTree())
    {
        std::cerr << "Unable to open output data file " << m_outputFileName << std::endl;
        return false;
    }

    return true;
}

bool DetectorApp::LoadConfig(const std::string& filename)
{
    std::cout<<"----------Detector Efficiency Calculation----------"<<std::endl;
    std::ifstream input(filename);
    if(!input.is_open())
    {
        std::cerr<<"Unable to open input config file "<<filename<<std::endl;
        return false;
    }

    std::string junk;
    std::string typeName;

    input >> junk >> m_inputFileName;
    input >> junk >> m_outputFileName;
    input >> junk >> m_deadChannelFileName;

    input >> junk >> m_nthreads;
    input >> junk >> typeName;

    std::cout << "Creating " << m_nthreads << " detector arrays..." << std::endl;
    ArrayType type = StringToArrayType(typeName);
    for(uint64_t i=0; i<m_nthreads; i++)
    {
        m_detectorList.push_back(CreateDetectorArray(type));
        if(m_deadChannelFileName != "None")
            m_detectorList.back()->SetDeadChannelMap(m_deadChannelFileName);
    }
    std::cout << "Done" << std::endl;

    std::cout << "Allocating " << m_nthreads << " threads..." << std::endl;
    m_resources = std::make_unique<Mask::ThreadPool<DetectorArray*, uint64_t>>(m_nthreads);
    std::cout << "Done" << std::endl;

    std::cout << "Opening input data file " << m_inputFileName << "..." << std::endl;
    m_fileReader.Open(m_inputFileName, "SimTree");
    if(!m_fileReader.IsOpen() || !m_fileReader.IsTree())
    {
        std::cerr << "Unable to open input data file " << m_inputFileName << std::endl;
        return false;
    }
    m_nentries = m_fileReader.GetSize();
    std::cout << "Done. Detected " << m_nentries << " events in the file." << std::endl;

    //Little bit of integer division mangling to make sure we read every event in file
    uint64_t quotient = m_nentries / m_nthreads;
    uint64_t remainder = m_nentries % m_nthreads;
    m_chunkSamples.push_back(quotient + remainder);
    for(uint64_t i=1; i<m_nthreads; i++)
        m_chunkSamples.push_back(quotient);

    std::cout << "Opening output data file " << m_outputFileName << std::endl;
    m_fileWriter.Open(m_outputFileName, "SimTree");
    if(!m_fileWriter.IsOpen() || !m_fileWriter.IsTree())
    {
        std::cerr << "Unable to open output data file " << m_outputFileName << std::endl;
        return false;
    }
    std::cout << "Done" << std::endl;

    std::cout << "Ready to launch!" << std::endl;

    return true;
}

void DetectorApp::Run()
{
	std::cout<<"Running efficiency calculation..."<<std::endl;

    if(m_detectorList.size() != m_nthreads)
    {
        std::cerr << "Detector list not equal to number of threads" << std::endl;
        return;
    }

    for(uint64_t i=0; i<m_detectorList.size(); i++)
    {
        //Create a job for the thread pool, using a lambda and providing a tuple of the arguments
        m_resources->PushJob({[this](DetectorArray* array, uint64_t chunkSamples)
            {
                if(system == nullptr)
		    	    return;

                std::vector<Mask::Nucleus> data;
                DetectorResult result;
	            for(uint64_t i=0; i<chunkSamples; i++)
	            {
                    m_fileReader.Read(data);
                    for(auto& nucleus : data)
                    {
                        result = array->IsDetected(nucleus);
                        if(result.detectFlag)
                        {
                            nucleus.isDetected = true;
                            nucleus.detectedKE = result.energy_deposited;
                            nucleus.detectedTheta = result.direction.Theta();
                            nucleus.detectedPhi = result.direction.Phi();
                            nucleus.detectedPos = result.direction;
                        }
                    }
		            m_fileWriter.PushData(data);
	            }
            },
        {m_detectorList[i], m_chunkSamples[i]} //arguments to function, in order
        }
        );
    }

	uint64_t size = m_fileReader.GetSize();
	uint64_t count = 0;
	double percent = 0.05;
	uint64_t flushVal = size*percent;
	uint64_t flushCount = 0;

	while(true)
	{
        if(count == flushVal)
	    {
	    	count = 0;
	    	++flushCount;
	    	std::cout<<"\rPercent of data written to disk: "<<percent*flushCount*100<<"%"<<std::flush;
	    }

		if(m_resources->IsFinished() && m_fileWriter.GetQueueSize() == 0)
				break;
		else if(m_fileWriter.Write())
				++count;
	}

    std::cout << std::endl;
	
}