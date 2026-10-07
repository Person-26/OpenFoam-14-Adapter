#ifndef FSI_PROPELLERLOAD_H
#define FSI_PROPELLERLOAD_H

#include "CouplingDataUser.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes the integrated actuator-disk propeller loads
//  (thrust force and reaction torque) from the OpenFOAM propellerDisk
//  fvModel to preCICE. One value is written per propeller hub, therefore
//  this data is meant to live on an interface with location type
//  "fixedPoints" registered at the propeller disk centres.
//
//  Sign convention: propellerDisk's normal is the direction of the thrust
//  on the propeller (it pushes the fluid the other way), and force() and
//  moment() are the fluid's force and moment on the propeller, so
//    - Thrust     = rho propellerDisk.force()
//    - PropTorque = rho propellerDisk.moment()
//  with rho the FSI density for an incompressible solver, whose disk loads
//  are kinematic (per unit density), and 1 for a compressible one.
class PropellerLoad : public CouplingDataUser
{
private:

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    //- Names of the propellerDisk fvModels to be coupled with
    //  (must match the order of the fixedPoints of the interface)
    std::vector<std::string> propellerNames_;

    //- Coupled quantity (true = moment, false = force)
    const bool isMoment_;

    //- Density the disk loads are scaled by (see above)
    const Foam::scalar rho_;

public:
    //- Constructor
    PropellerLoad(
        const Foam::fvMesh& mesh,
        std::vector<std::string> propellerNames,
        bool isMoment,
        Foam::scalar rho = 1.0);

    //- Write the propeller loads to the buffer
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the propeller loads from the buffer (not implemented)
    void read(double* buffer, const unsigned int dim) final;

    //- For this data set, only the fixedPoints location type is supported
    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif