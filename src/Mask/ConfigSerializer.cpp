#include "ConfigSerializer.h"
#include "yaml-cpp/yaml.h"
#include "RxnType.h"

namespace Mask {

    static void SerializeDecay(YAML::Emitter& yamlStream, const StepParameters& params)
    {
        yamlStream << YAML::BeginMap;
		yamlStream << YAML::Key << "Type" << RxnTypeToString(params.rxnType); 
		yamlStream << YAML::Key << "Reactants" << YAML::Value << YAML::BeginSeq;
		for (std::size_t i=0; i<params.Z.size(); i++)
		{
			yamlStream << YAML::BeginMap;
			yamlStream << YAML::Key << "Z" << YAML::Value << params.Z[i];
			yamlStream << YAML::Key << "A" << YAML::Value << params.A[i];
			yamlStream << YAML::EndMap;
		}
		yamlStream << YAML::EndSeq;
        yamlStream << YAML::Key << "PhiMin(deg)" << YAML::Value << params.phiMin;
        yamlStream << YAML::Key << "PhiMax(deg)" << YAML::Value << params.phiMax;
        yamlStream << YAML::Key << "ResidualExcitationMean(MeV)" << YAML::Value << params.meanResidualEx;
        yamlStream << YAML::Key << "ResidualExcitationSigma(MeV)" << YAML::Value << params.sigmaResidualEx;
		yamlStream << YAML::Key << "AngularDistributionFile" << YAML::Value << params.angularDistFile;
		yamlStream << YAML::EndMap;
    }

    static StepParameters DeserializeDecay(YAML::Node& yamlStream)
    {
        StepParameters params;
        params.rxnType = StringToRxnType(yamlStream["Type"].as<std::string>());
        YAML::Node nuclei = yamlStream["Reactants"];
        for(auto nucleus : nuclei)
        {
            params.Z.push_back(nucleus["Z"].as<int>());
            params.A.push_back(nucleus["A"].as<int>());
        }
        params.phiMin = yamlStream["PhiMin(deg)"].as<double>();
        params.phiMax = yamlStream["PhiMax(deg)"].as<double>();
        params.meanResidualEx = yamlStream["ResidualExcitationMean(MeV)"].as<double>();
        params.sigmaResidualEx = yamlStream["ResidualExcitationSigma(MeV)"].as<double>();
        params.angularDistFile = yamlStream["AngularDistributionFile"].as<std::string>();
        return params;
    }

    static void SerializeReaction(YAML::Emitter& yamlStream, const StepParameters& params)
    {
        yamlStream << YAML::BeginMap;
		yamlStream << YAML::Key << "Type" << RxnTypeToString(params.rxnType); 
		yamlStream << YAML::Key << "Reactants" << YAML::Value << YAML::BeginSeq;
		for (std::size_t i=0; i<params.Z.size(); i++)
		{
			yamlStream << YAML::BeginMap;
			yamlStream << YAML::Key << "Z" << YAML::Value << params.Z[i];
			yamlStream << YAML::Key << "A" << YAML::Value << params.A[i];
			yamlStream << YAML::EndMap;
		}
		yamlStream << YAML::EndSeq;
        yamlStream << YAML::Key << "BeamEnergyMean(MeV)" << YAML::Value << params.meanBeamEnergy;
        yamlStream << YAML::Key << "BeamEnergySigma(MeV)" << YAML::Value << params.sigmaBeamEnergy;
        yamlStream << YAML::Key << "ThetaType" << YAML::Value << RxnThetaTypeToString(params.thetaType);
        yamlStream << YAML::Key << "ThetaMin(deg)" << YAML::Value << params.thetaMin;
        yamlStream << YAML::Key << "ThetaMax(deg)" << YAML::Value << params.thetaMax;
        yamlStream << YAML::Key << "PhiMin(deg)" << YAML::Value << params.phiMin;
        yamlStream << YAML::Key << "PhiMax(deg)" << YAML::Value << params.phiMax;
        yamlStream << YAML::Key << "ResidualExcitationMean(MeV)" << YAML::Value << params.meanResidualEx;
        yamlStream << YAML::Key << "ResidualExcitationSigma(MeV)" << YAML::Value << params.sigmaResidualEx;
		yamlStream << YAML::EndMap;
    }

    static StepParameters DeserializeReaction(YAML::Node& yamlStream)
    {
        StepParameters params;
        params.rxnType = StringToRxnType(yamlStream["Type"].as<std::string>());
        YAML::Node nuclei = yamlStream["Reactants"];
        for(auto nucleus : nuclei)
        {
            params.Z.push_back(nucleus["Z"].as<int>());
            params.A.push_back(nucleus["A"].as<int>());
        }
        params.meanBeamEnergy = yamlStream["BeamEnergyMean(MeV)"].as<double>();
        params.sigmaBeamEnergy = yamlStream["BeamEnergySigma(MeV)"].as<double>();
        params.thetaType = StringToRxnThetaType(yamlStream["ThetaType"].as<std::string>());
        params.thetaMin = yamlStream["ThetaMin(deg)"].as<double>();
        params.thetaMax = yamlStream["ThetaMax(deg)"].as<double>();
        params.phiMin = yamlStream["PhiMin(deg)"].as<double>();
        params.phiMax = yamlStream["PhiMax(deg)"].as<double>();
        params.meanResidualEx = yamlStream["ResidualExcitationMean(MeV)"].as<double>();
        params.sigmaResidualEx = yamlStream["ResidualExcitationSigma(MeV)"].as<double>();
        return params;
    }

    static void SerializeTarget(YAML::Emitter& yamlStream, const LayeredTarget& target)
    {

        yamlStream << YAML::BeginSeq;
        for (int i=0; i<target.GetNumberOfLayers(); i++)
        {
            const Target& layer = target.GetLayerInfo(i);
            yamlStream << YAML::BeginMap;
            yamlStream << YAML::Key << "Thickness(ug/cm^2)" << YAML::Value << layer.GetThickness();
            yamlStream << YAML::Key << "Elements" << YAML::Value << YAML::BeginSeq;
            for (int j=0; j<layer.GetNumberOfElements(); j++)
            {
                yamlStream << YAML::BeginMap;
                yamlStream << YAML::Key << "Z" << YAML::Value << layer.GetElementZ(j);
                yamlStream << YAML::Key << "A" << YAML::Value << layer.GetElementA(j);
                yamlStream << YAML::Key << "S" << YAML::Value << layer.GetElementStoich(j);
                yamlStream << YAML::EndMap;
            }
            yamlStream << YAML::EndSeq;
            yamlStream << YAML::EndMap;
        }
        yamlStream << YAML::EndSeq;
    }

    static LayeredTarget DeserializeTarget(YAML::Node& yamlStream)
    {
        LayeredTarget target;
        std::vector<int> z, a, s;
        double thickness;

        for (auto layer : yamlStream)
        {
            thickness = layer["Thickness(ug/cm^2)"].as<double>();
            YAML::Node elements = layer["Elements"];
            z.clear();
            a.clear();
            s.clear();
            for(auto element : elements)
            {
                z.push_back(element["Z"].as<int>());
                a.push_back(element["A"].as<int>());
                s.push_back(element["S"].as<int>());
            }
            target.AddLayer(z, a, s, thickness);
        }

        return target;
    }

    bool ConfigSerializer::SerializeConfig(const std::string& configfile, const AppParameters& params)
    {
		std::ofstream output(configfile);
		if (!output.is_open())
		{
			std::cerr << "Could not open output file " << configfile << std::endl;
			return false;
		}

		YAML::Emitter yamlStream;
		yamlStream << YAML::BeginMap;
		yamlStream << YAML::Key << "OutputFile" << YAML::Value << params.outputFileName;
		yamlStream << YAML::Key << "Threads" << YAML::Value << params.nThreads;
		yamlStream << YAML::Key << "ReactionSamples" << YAML::Value << params.nSamples;
		yamlStream << YAML::Key << "ReactionChain" << YAML::Value << YAML::BeginSeq;
		for (auto& step : params.chainParams)
		{
			switch(step.rxnType)
            {
                case RxnType::Reaction: SerializeReaction(yamlStream, step); break;
                case RxnType::Decay: SerializeDecay(yamlStream, step); break;
                case RxnType::None:
                {
                    std::cerr << "Error serializing config: none type reaction found!" << std::endl;
                    return false;
                }
            }
		}
        yamlStream << YAML::EndSeq;
        yamlStream << YAML::Key << "TargetLayers" << YAML::Value;
        SerializeTarget(yamlStream, params.target);

        output << yamlStream.c_str();
        output.close();
        return true;
    }

    bool ConfigSerializer::DeserializeConfig(const std::string& configfile, AppParameters& params)
    {
        YAML::Node data;
		try 
		{
			data = YAML::LoadFile(configfile);
		}
		catch(YAML::ParserException& e)
		{
			std::cerr << "Error reading config " << configfile << ": " << e.what() <<std::endl;
			return false;
		}

        params.outputFileName = data["OutputFile"].as<std::string>();
        params.nThreads = data["Threads"].as<uint32_t>();
        params.nSamples = data["ReactionSamples"].as<uint64_t>();

        auto steps = data["ReactionChain"];
        for (auto step : steps)
        {
            RxnType type = StringToRxnType(step["Type"].as<std::string>());
            switch(type)
            {
                case RxnType::Decay:
                {
                    params.chainParams.push_back(DeserializeDecay(step));
                    break;
                }
                case RxnType::Reaction:
                {
                    params.chainParams.push_back(DeserializeReaction(step));
                    break;
                }
                case RxnType::None:
                {
                    std::cerr << "Error deserializing config: None-type reaction found" << std::endl;
                    return false;
                }
            }
        }
        auto layers = data["TargetLayers"];
        params.target = DeserializeTarget(layers);
        return true;
    }
}