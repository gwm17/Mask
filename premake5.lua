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

	filter "system:windows"
		systemversion "latest"
		--ROOT cannot be located using the config tools, so we must query for the ROOTSYS env variable
		--Have to use an if statement to hide this (@penguin for example doesn't have a ROOTSYS)
		if os.host() == windows then
			rootpath = os.getenv("ROOTSYS")

			includedirs {
				"include",
				"./",
				rootpath .. "include"
			}

			links {
				rootpath .. "lib/**.lib"
			}
		end

	filter "system:macosx or linux"
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
		"include/*.h"
	}

	includedirs {
		"include"
	}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"