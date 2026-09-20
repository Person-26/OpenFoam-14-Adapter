#ifndef FSI_PROPELLERRPM_H
#define FSI_PROPELLERRPM_H

#include "CouplingDataUser.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that reads a commanded rotation speed (scalar, one per propeller
//  hub) from preCICE and applies it to the OpenFOAM propellerDisk fvModels.
//  The sign of the RPM selects the rotation direction, which lets a coupled
//  controller command both the magnitude and the sense (e.g. alternating
//  spin directions on a multirotor).
//
//  This data is meant to live on an interface with location type
//  "fixedPoints" whose vertices coincide with the propeller disk centres
//  (one vertex per propeller, in the same order as the propellers list).
//
//  Requires the propellerDisk fvModel to expose setRotationSpeed() (a small
//  OpenFOAM-14 extension); the disk's optional "adjustment" PI controller
//  must be disabled for external control to take effect.
class PropellerRPM : public CouplingDataUser
{
private:

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    //- Names of the propellerDisk fvModels to be coupled with
    //  (must match the order of the fixedPoints of the interface)
    std::vector<std::string> propellerNames_;

public:
    //- Constructor
    PropellerRPM(
        const Foam::fvMesh& mesh,
        std::vector<std::string> propellerNames);

    //- Write the RPM to the buffer (not implemented)
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the commanded RPM from the buffer and set the disk rotation speed
    void read(double* buffer, const unsigned int dim) final;

    //- For this data set, only the fixedPoints location type is supported
    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
