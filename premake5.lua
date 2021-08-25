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
		"src/**.cpp",
		"include/**.h",
		"src/**.cxx"
	}

	includedirs {
		"include",
		"./"
	}

	buildoptions {
		"`root-config --cflags`"
	}

	linkoptions {
		"`root-config --glibs`"
	}

	prebuildcommands {
		"rootcint -f src/kinematics_dict.cxx include/Kinematics.h include/LinkDef_Kinematics.h",
		"{MOVE} src/kinematics_dict_rdict.pcm bin/"
	}