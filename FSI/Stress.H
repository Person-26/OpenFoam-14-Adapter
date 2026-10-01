#ifndef FSI_STRESS_H
#define FSI_STRESS_H

#include "ForceBase.H"

#include <utility>
#include <vector>

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads stress [N/m^2]:
// This is essentially a force (in spatial coordinates) scaled by the deformed
// cell face. Thus, a consistent quantity. Calculation concept has been copied
// from the force module, but the scaled version here is commonly used in FEM
// applications.
class Stress : public ForceBase
{

private:
    //- Stress field
    Foam::volVectorField* Stress_;

    //- Buffer locations of the two sides of each zero-thickness baffle face
    //  in the interface (coincident face centres, opposing normals)
    std::vector<std::pair<std::size_t, std::size_t>> bafflePairs_;

    //- Whether bafflePairs_ has been built
    bool bafflePairsBuilt_ = false;

    //- Find the baffle face pairs among the interface faces
    void findBafflePairs();

public:
    //- Constructor
    Stress(
        const Foam::fvMesh& mesh,
        const std::string solverType);

    //- Write the stress values into the buffer
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the stress values from the buffer
    void read(double* buffer, const unsigned int dim) final;

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;

    //- Returns the face normal vectors (no multiplication by area)
    Foam::tmp<Foam::vectorField> getFaceVectors(const unsigned int patchID) const final;

    //- Destructor
    ~Stress();
};
}
}

#endif
