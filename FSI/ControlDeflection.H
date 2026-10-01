#ifndef FSI_CONTROLDEFLECTION_H
#define FSI_CONTROLDEFLECTION_H

#include "CouplingDataUser.H"

namespace preciceAdapter
{
namespace FSI
{

//- Configuration of a single (rigidly rotating) control surface
struct ControlSurfaceConfig
{
    std::string patch;          // OpenFOAM patch of the surface
    std::vector<double> hinge;  // hinge point on the rotation axis [m]
    std::vector<double> axis;   // unit rotation axis [-]

    //- Last commanded deflection angle received for this surface [rad].
    //  Shared between the ControlDeflection reader and the ControlCommand
    //  writer so the latter can relay the command on to the structural solver.
    double lastDeflection = 0.0;

    //- Set when a ServoAngle reader is configured: the patch then follows the
    //  angle the structural servo reached, and the Deflection reader only
    //  records the command.
    bool servoDriven = false;
};

//- Class that reads a per-surface deflection angle (scalar) from preCICE and
//  rotates the corresponding control-surface patch about its hinge. It reads
//  either the commanded angle ("Deflection") or the angle an actuated servo
//  reached ("ServoAngle"); when both are read, the servo angle moves the patch
//  and the command is only recorded. The
//  resulting displacement is written into the pointDisplacement field on that
//  patch, from which OpenFOAM's dynamicMotionSolverFvMesh (with
//  displacementLaplacian) propagates the deformation into the fluid mesh.
//
//  This data is meant to live on an interface with location type
//  "fixedPoints" whose vertices coincide with the surface hinges (one vertex
//  per control surface, in the same order as the controlSurfaces list).
class ControlDeflection : public CouplingDataUser
{
private:

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    //- Configuration of the control surfaces (patch + hinge + axis). Owned by
    //  the FSI object and shared with the ControlCommand writer.
    std::vector<ControlSurfaceConfig>* controlSurfaces_;

    //- pointDisplacement field name
    std::string namePointDisplacement_;

    //- Whether this reader takes the servo angle instead of the command
    bool servoAngle_;

    //- Reference (initial) points of each control-surface patch, stored at
    //  initialize() so the rotation is always measured from the reference
    //  (undeflected) configuration.
    std::vector<std::vector<double>> referencePoints_;

    //- A hinge-line point of a control surface that it shares with another
    //  moving patch (e.g. a deforming plate the surface is hinged to)
    struct HingeLink
    {
        double span;          // position along the hinge axis [m]
        Foam::label patchID;  // the other patch
        Foam::label index;    // the point's index in that patch
    };

    //- Hinge links of each control surface, sorted by span. The surface is
    //  translated with the displacement of these points, so it stays attached
    //  along the hinge when the patch it is hinged to deforms.
    std::vector<std::vector<HingeLink>> hingeLinks_;

    //- Initialize the reference points of all control-surface patches
    void storeReferencePoints();

    //- Find the hinge-line points each control surface shares with other
    //  patches
    void findHingeLinks();

public:
    //- Constructor
    ControlDeflection(
        const Foam::fvMesh& mesh,
        std::vector<ControlSurfaceConfig>& controlSurfaces,
        bool servoAngle = false);

    //- Write the deflection to the buffer (not implemented)
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the deflection angles from the buffer and rotate the patches
    void read(double* buffer, const unsigned int dim) final;

    //- For this data set, only the fixedPoints location type is supported
    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;

    //- Initialize the reference points (called after configuration)
    void initialize() final;
};

}
}

#endif