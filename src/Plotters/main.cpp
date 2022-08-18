#include "RootPlotter.h"
#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        std::cerr<<"Root plotter requires two commandline arguments: path to input file and path to outputfile"<<std::endl;
        return 1;
    }

    RootPlotter plotter;
    plotter.Run(argv[1], argv[2]);
}