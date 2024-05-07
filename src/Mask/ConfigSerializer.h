#ifndef CONFIG_SERIALIZER_H
#define CONFIG_SERIALIZER_H

#include "MaskApp.h"

namespace Mask {

    class ConfigSerializer
    {
    public:
        static bool SerializeConfig(const std::string& configfile, const AppParameters& params);
        static bool DeserializeConfig(const std::string& configfile, AppParameters& params);
    };
}

#endif