# MASK: Monte cArlo Simulation of Kinematics
MASK is a Monte Carlo simulation of reaction kinematics for use detector systems at Florida State University.
MASK is capable of simulating multi-step kinematic reaction-decay sequences, storing data in a lightweight binary format, after which the kinematic data can be fed to a detector geometry for efficiency testing. Currently geometries for ANASEN and SABRE are included in the code.

## Building MASK
Dowload the repository from github. The code is to be built using Premake5, an open source and free project building software. To build on Linux and MacOSX run
`premake5 gmake2`
to generate the makefiles. Then execute
`make`
to build the program. To build on Windows run
`premake5 vs20xx`
where xx should be replaced with your version of Visual Studio, and then build as a normal Visual Studio project. For more documentation see the Premake wiki. By default `make` runs in Debug mode. For release mode use `make config=release`

### Building RootPlot
One of the executables, RootPlot, requires linking against external libraries from the ROOT cern analysis package. Path to the necessary header files and libraries must be set by the user on a machine-by-machine basis in the premake5.lua file. Set the ROOTIncludepath and ROOTLibpath to match your install. An easy way to check for the paths on a unix system is through the root-config tool. root-config --cflags has the include path (the path after -I) and root-config --glibs has the lib path (after -L). If you do not have ROOT installed, you will not be able to compile RootPlot, and you will not be able to use the generic `make` command, which will try to build all executables. Instead use `make Mask` or `make DetectEff` (or `make config=release Mask`, etc).

## Running MASK
By default MASK is capable of simulating reactions of up to three steps. Here is a brief outline of each type:

0. A reaction of type 0 is a pure decay. It is assumed isotropic by default; any other case will require the modification of the code.
1. A reaction of type 1 is a pure reaction. It can incorporate all of the input file sampling parameters.
2. A reaction of type 2 is a reaction followed by a subsequent decay of the residual nucleus. Again, all sampling is allowed.
3. A reaction of type 3 is a reaction followed by a subsequent decay of the residual, followed by a decay of one of the products. Again, all sampling is allowed

For decays, a specific angular distribution can be given as input as a text file with values of coefficiencts of a Legendre polynomial series. Examples can be found in the `./etc` directory, including an isotropic case. It is assumed that the decays in the center-of-mass frame are isotropic in phi (i.e. m=0). Decay1 corresponds to the first decay, if there are multiple steps, Decay2 to the second. If there are no decays, these parameters are not used (or if only one decay, Decay2_AngularMomentum is not used). The input file requires that the user include target information, which will be used to calculate energy loss for all of the reactants and reaction products. The target can contain layers, and each layer can be composed of a compound of elements with a given stoichiometry. If the user wishes to not include energy loss in the kinematics, simply give all target layers a thickness of 0. Note that more layers and more thickness = more time spent calculating energy loss. These energy loss methods are only applicable for solid targets, and should not be applied to gas or liquid targets. Energy loss calculations have a stated uncertainty of approximately five percent.

To choose which detector scheme is run, modify the main function in DetectorEfficiency.cpp. The included geometries also have options to do an internal geometry consistency check and print out coordinates for drawing the detector arrays.

To run MASK simply do the following from the MASK directory:

`./bin/Mask input.txt`

Input.txt can be replaced by any text file with the correct format.

To run DetEff use the format

`./bin/DetEff <kinematics_datafile> <new_detection_datafile> <new_detection_statsfile>`

where the detection datafile contains all of the kinematics data as well as information about which particles are detected (this is in the mask file format) and the statsfile is a text file containing efficiency statistics.

RootPlot is run as

`./bin/RootPlot <datafile> <outputfile>`

where the datafile can be either the datafile from Mask or the datafile from DetEff. The outputfile is saved in the ROOT file format.

PyPlotter is run as

`./bin/PyPlotter <datafile> <outputfile>`

where again datafile can come from either Mask or DetEff. The outputfile is saved in the Python pickle file format. Pickled files can be reopened using the PyPlotViewer script run as
`./bin/PyPlotViewer <picklefile>`

## Requirements
MASK, the kinematics simulation, and DetEff, the detection efficiency simulation, require no external dependancies. The RootPlot plotting tool requires the ROOT cern data analysis package, and the PyPlotter tool requires Python3 and the matplotlib and numpy packages. 