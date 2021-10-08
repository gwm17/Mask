workspace "Mask" 
	configurations {
		"Debug",
		"Release"
	}

project "Mask"
	kind "ConsoleApp"
	language "C++"
	targetdir "bin"
	objdir "objs"
	cppdialect "C++11"

	files {
		"src/*.cpp",
		"include/*.h"
	}

	includedirs {
		"include"
	}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"

project	"RootPlot"
	kind "ConsoleApp"
	language "C++"
	targetdir "bin"
	objdir "objs"
	cppdialect "c++11"

	files {
		"src/Plotters/ROOT/RootPlotter.cpp",
		"src/MaskFile.cpp",
		"src/Nucleus.cpp",
		"src/Vec4.cpp",
		"src/Vec3.cpp",
		"src/MassLookup.cpp",
		"include/*.h"
	}

	--User specified path to ROOT CERN libraries--
	ROOTIncludepath = "/home/gordon/cern/root-6.22.02/root-install/include"
	ROOTLibpath = "/home/gordon/cern/root-6.22.02/root-install/lib"

	includedirs {
		"include"
	}

	sysincludedirs {
		ROOTIncludepath
	}

	libdirs {
		ROOTLibpath
	}

	links {
		"Gui", "Core", "Imt", "RIO", "Net", "Hist", 
		"Graf", "Graf3d", "Gpad", "ROOTDataFrame", "ROOTVecOps",
		"Tree", "TreePlayer", "Rint", "Postscript", "Matrix",
		"Physics", "MathCore", "Thread", "MultiProc", "m", "dl"
	}

	filter "system:macosx or linux"
		linkoptions {
			"-pthread"
		}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"

project	"DetectEff"
	kind "ConsoleApp"
	language "C++"
	targetdir "bin"
	objdir "objs"
	cppdialect "c++11"

	files {
		"src/Detectors/*.cpp",
		"src/MaskFile.cpp",
		"src/Nucleus.cpp",
		"src/Vec4.cpp",
		"src/Vec3.cpp",
		"src/MassLookup.cpp",
		"src/Rotation.cpp",
		"src/Target.cpp",
		"src/EnergyLoss.cpp",
		"src/RandomGenerator.cpp",
		"include/*.h"
	}

	includedirs {
		"include"
	}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"