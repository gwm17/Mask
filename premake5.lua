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

	prebuildcommands {
		"rootcint -f src/kinematics_dict.cxx include/Kinematics.h include/LinkDef_Kinematics.h",
		"{MOVE} src/kinematics_dict_rdict.pcm bin/"
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