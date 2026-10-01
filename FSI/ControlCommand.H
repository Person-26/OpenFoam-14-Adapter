#ifndef FSI_CONTROLCOMMAND_H
#define FSI_CONTROLCOMMAND_H

#include "ControlDeflection.H"
#include "CouplingDataUser.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes the commanded control-surface angles (scalar, one per
//  surface) by relaying the deflection that the ControlDeflection reader last
//  applied. It lets the fluid forward the flight-controller command to the
//  structural solver's servo hinge without a third coupling link between the
//  kinematics and the solid participants.
//
//  This data lives on the same "fixedPoints" interface as ControlDeflection,
//  with vertices at the surface hinges in the controlSurfaces order.
class ControlCommand : public CouplingDataUser
{
private:

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    //- Control-surface configuration, owned by the FSI object and shared with
    //  the ControlDeflection reader (which stores the applied angle in it).
    std::vector<ControlSurfaceConfig>* controlSurfaces_;

public:
    //- Constructor
    ControlCommand(
        const Foam::fvMesh& mesh,
        std::vector<ControlSurfaceConfig>& controlSurfaces);

    //- Write the last applied per-surface deflection to the buffer
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read is not implemented (this data is written to preCICE, not read)
    void read(double* buffer, const unsigned int dim) final;

    //- For this data set, only the fixedPoints location type is supported
    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;

    //- Initialize (no-op; the command is shared through the configuration)
    void initialize() final;
};

}
}

#endif
