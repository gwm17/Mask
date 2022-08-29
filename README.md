# MASK: Monte cArlo Simulation of Kinematics
MASK is a Monte Carlo simulation of reaction kinematics for use detector systems at Florida State University.
MASK is capable of simulating multi-step kinematic reaction-decay sequences, storing data in ROOT Trees, after which the kinematic data can be fed to a detector geometry for efficiency testing. Currently geometries for ANASEN and SABRE are included in the code.

## Building MASK
Dowload the repository from github. CMake is use to build the project; in most environments you can build Mask using the following methods:
- `mkdir build`
- `cd build`
- `cmake ..`
- `make`

Executables will be installed to the repostiory's `bin` directory. Libraries will be installed to the `lib` directory.

By default Mask builds for release. To build for debug replace `cmake ..` with `cmake -DCMAKE_BUILD_TYPE=Debug ..`. Mask uses CMake to find the installed ROOT libraries and headers.

## Using the kinematics simulation
By default MASK is capable of simulating reactions of up to three steps. Here is a brief outline of each type:

0. A pure decay involves a "target" decaying into an ejectile and residual. It is assumed isotropic by default; any other case will require the modification of the code.
1. A single step reaction involves a target being hit by a projectile and emitting an ejectile and a residual. It can incorporate all of the input file sampling parameters.
2. A two step reaction is a single step reaction where the residual subsequently decays into two daughters (called breakups). Again, all sampling is allowed.
3. A three step reaction is a two step reaction where one of the breakup particles subsequently decays into two more breakup particles. Again, all sampling is allowed

For decays, a specific angular distribution can be given as input as a text file with values of coefficiencts of a Legendre polynomial series. Examples can be found in the `./etc` directory, including an isotropic case. It is assumed that the decays in the center-of-mass frame are isotropic in phi (i.e. m=0). Decay1 corresponds to the first decay, if there are multiple steps, Decay2 to the second. If there are no decays, these parameters are not used (or if only one decay, Decay2_AngularMomentum is not used). The input file requires that the user include target information, which will be used to calculate energy loss for all of the reactants and reaction products. The energy loss through materials is calculated using the `catima` library (found in `src/vendor`), which is a C/C++ interface to the `atima` library (the same energy loss methods used by LISE). The target can contain layers, and each layer can be composed of a compound of elements with a given stoichiometry. If the user wishes to not include energy loss in the kinematics, simply give all target layers a thickness of 0. Note that more layers and more thickness = more time spent calculating energy loss. These energy loss methods are only applicable for solid targets, and should not be applied to gas or liquid targets. Energy loss calculations have a stated uncertainty of approximately five percent.

The input file requires that the user include target information, which will be used to calculate energy loss for all of the reactants and reaction products. The target can contain layers, and each layer can be composed of a compound of elements with a given stoichiometry. If the user wishes to not include energy loss in the kinematics, simply give all target layers a thickness of 0. Note that more layers and more thickness = more time spent calculating energy loss. These energy loss methods are only applicable for solid targets, and should not be applied to gas or liquid targets. Energy loss calculations have a stated uncertainty of approximately five percent. The energy loss library used is called `catima` and can be found [here](https://github.com/gwm17/catima).

To run MASK simply do the following from the MASK repository:

`./bin/MaskApp input.txt`

Input.txt can be replaced by any text file with the correct format.

## Using the detector geometry simulation
Detector geometry is encoded using ROOT math libraries in the `src/Detectors` folder. Two different detector geometries are already present: SPS-SABRE and ANASEN. To add a new geometry, follow the guidelines outlined by each of these cases.

To choose which detector scheme is run, modify the main function in `src/Detectors/main.cpp`. The included geometries also have options to do an internal geometry consistency check and print out coordinates for drawing the detector arrays, which can be useful for testing.



To run DetEff use the format

`./bin/DetEff <kinematics_datafile> <new_detection_datafile> <new_detection_statsfile>`

where the detection datafile contains all of the kinematics data as well as information about which particles are detected and the statsfile is a text file containing efficiency statistics.

## Data visualization
All data is saved as ROOT trees of std::vectors of Mask::Nucleus classes. To enable this, a ROOT dictionary is generated and linked into a shared library found in the `lib` directory of the repository. This allows the user to link to the shared library for accessing and analyzing the data generated by MASK.

Mask also provides a default visualization tool called RootPlot. RootPlot is run as

`./bin/RootPlot <datafile> <outputfile>`

where the datafile can be either the datafile from Mask or the datafile from DetEff. The outputfile is saved in the ROOT file format.

## Requirements
ROOT version 6.22 or greater is required
CMake version 3.0 or greater is required

