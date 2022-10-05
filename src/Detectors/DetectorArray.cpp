#include "DetectorArray.h"
#include "AnasenArray.h"
#include "SabreArray.h"

DetectorArray* CreateDetectorArray(ArrayType type)
{
    switch(type)
    {
        case ArrayType::None: return nullptr;
        case ArrayType::Anasen: return new AnasenArray();
        case ArrayType::Sabre: return new SabreArray();
    }
    return nullptr;
}

std::string ArrayTypeToString(ArrayType type)
{
    switch(type)
    {
        case ArrayType::None: return "None";
        case ArrayType::Anasen: return "Anasen";
        case ArrayType::Sabre: return "Sabre";
    }
    return "None";
}

ArrayType StringToArrayType(const std::string& value)
{
    if (value == "Anasen")
        return ArrayType::Anasen;
    else if (value == "Sabre")
        return ArrayType::Sabre;
    else
        return ArrayType::None;
}