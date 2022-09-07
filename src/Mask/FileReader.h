#ifndef FILE_READER_H
#define FILE_READER_H

#include "Nucleus.h"

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <string>
#include <mutex>
#include <atomic>

namespace Mask {

    class FileReader
    {
    public:
        FileReader();
        FileReader(const std::string& filename, const std::string& treename);
        ~FileReader();

        void Open(const std::string& filename, const std::string& treename); //Not thread safe
        void Close(); //Not thread safe

        /*
            Read: fills entry to given dataHandle. Returns true if data was successfully filled, otherwise returns false
        */
        bool Read(std::vector<Nucleus>& dataHandle); //Thread safe
        uint64_t GetSize() { return m_size; }//In entries (implicitly thread safe)
        bool IsOpen() { return m_file == nullptr ? false : m_file->IsOpen(); } //Should be safe?
        bool IsTree() { return m_tree != nullptr; } //Should be safe?

    private:
        TFile* m_file;
        TTree* m_tree;

        std::vector<Nucleus>* m_branchHandle;

        std::mutex m_fileMutex;
        std::atomic<uint64_t> m_currentEntry;
        std::atomic<uint64_t> m_size; //in entries
    };
}

#endif