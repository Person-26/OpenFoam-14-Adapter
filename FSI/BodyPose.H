#ifndef FSI_BODYPOSE_H
#define FSI_BODYPOSE_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FSI
{

//- Reads the vehicle's pose from preCICE: its translation
//  ("BodyDisplacement") or rotation vector ("BodyRotation") relative to the
//  start, one vector on a single fixed point. It updates the pose that the
//  bodyDisplacementLaplacian mesh mover, the propeller disks and the load
//  writers use (see BodyFrame.H).
class BodyPose : public CouplingDataUser
{
private:
    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    //- Registry field this reader sets ("bodyDisplacement" or "bodyRotation")
    Foam::word fieldName_;

    //- The latest four samples, oldest first (a sample arrives once per
    //  coupling window; see read() for the interpolation), the latest's
    //  arrival time and the window length
    Foam::vector samples_[4] = {Foam::Zero, Foam::Zero, Foam::Zero, Foam::Zero};
    Foam::scalar tLatest_ = 0;
    Foam::scalar window_ = 0;
    bool haveSample_ = false;

    //- Patches of the vehicle's surface that take its rigid-body velocity
    //  (FSI dictionary: bodyWallPatches; set by the rotation reader, which
    //  reads after the displacement one)
    Foam::wordList wallPatches_;

    //- Set the rigid-body velocity on the wall patches
    void setWallVelocity() const;

public:
    //- Constructor
    BodyPose(const Foam::fvMesh& mesh, bool rotation);

    //- Write (not implemented)
    std::size_t write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the pose vector
    void read(double* buffer, const unsigned int dim) final;

    //- A single fixed point
    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;
};

}
}

#endif
