workspace "Mask" 
	configurations {
		"Release",
		"Debug"
	}

project "Mask"
	kind "StaticLib"
	language "C++"
	targetdir "lib"
	objdir "objs"
	cppdialect "C++17"

	files {
		"src/Mask/*.cpp",
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
	cppdialect "c++17"

	files {
		"src/Plotters/ROOT/RootPlotter.cpp",
		"include/*.h"
	}

	--User specified path to ROOT CERN libraries--
	ROOTIncludepath = "/Users/gordon/Cern/root/include"
	ROOTLibpath = "/Users/gordon/Cern/root/lib"

	includedirs {
		"include"
	}

	sysincludedirs {
		ROOTIncludepath
	}

	libdirs {
		"lib/",
		ROOTLibpath
	}

	links {
		"Mask", "Gui", "Core", "Imt", "RIO", "Net", "Hist", 
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

project "MaskApp"
	kind "ConsoleApp"
	language "C++"
	targetdir "bin"
	objdir "objs"
	cppdialect "c++17"

	files {
		"src/MaskApp/*.cpp"
	}

	includedirs {
		"include"
	}

	links {
		"Mask"
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
	cppdialect "c++17"

	files {
		"src/Detectors/*.cpp",
		"include/*.h"
	}

	includedirs {
		"include"
	}

	links {
		"Mask"
	}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"
