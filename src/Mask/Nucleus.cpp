#include "Nucleus.h"

namespace Mask {

    Nucleus CreateNucleus(int z, int a)
    {
        Nucleus nuc;
        nuc.Z = z;
        nuc.A = a;
        nuc.groundStateMass = MassLookup::GetInstance().FindMass(z, a);
		nuc.isotopicSymbol = MassLookup::GetInstance().FindSymbol(z, a);
		nuc.vec4 = ROOT::Math::PxPyPzEVector(0., 0., 0., nuc.groundStateMass);
        return nuc;
    }

    bool EnforceDictionaryLinked() { return true; }
}