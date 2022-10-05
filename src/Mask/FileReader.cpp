#include "FileReader.h"

namespace Mask {

    FileReader::FileReader() :
        m_file(nullptr), m_tree(nullptr), m_branchHandle(nullptr), m_currentEntry(0), m_size(0)
    {
    }

    FileReader::FileReader(const std::string& filename, const std::string& treename) :
        m_file(nullptr), m_tree(nullptr), m_branchHandle(nullptr), m_currentEntry(0), m_size(0)
    {
        Open(filename, treename);
    }

    FileReader::~FileReader()
    {
        Close();
    }

    void FileReader::Open(const std::string& filename, const std::string& treename)
    {
        if(m_file != nullptr || m_tree != nullptr)
            Close();

        m_file = TFile::Open(filename.c_str(), "READ");
        if(m_file != nullptr && m_file->IsOpen())
        {
            m_tree = (TTree*) m_file->Get(treename.c_str());
            if(m_tree == nullptr)
            {
                m_file->Close();
                delete m_file;
                m_file = nullptr;
            }
            m_branchHandle = new std::vector<Nucleus>();
            m_tree->SetBranchAddress("nuclei", &m_branchHandle);
            m_size = m_tree->GetEntries();
            m_currentEntry = 0; //Reset file position
        }
    }

    void FileReader::Close()
    {
        if(m_file != nullptr && m_file->IsOpen())
        {
            m_file->Close();
            delete m_file;
            m_file = nullptr;
        }
    }

    bool FileReader::Read(std::vector<Nucleus>& dataHandle)
    {
        std::scoped_lock<std::mutex> guard(m_fileMutex);
        int bytes = m_tree->GetEntry(m_currentEntry);
        if(bytes != 0)
        {
            dataHandle = *m_branchHandle;
            m_currentEntry++;
            return true;
        }

        return false;
    }
}