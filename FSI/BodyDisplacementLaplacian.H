#ifndef FSI_BODYDISPLACEMENTLAPLACIAN_H
#define FSI_BODYDISPLACEMENTLAPLACIAN_H

#include "displacementLaplacian_fvMotionSolver.H"

namespace Foam
{
namespace fvMotionSolvers
{

//- Mesh motion of a vehicle's fluid domain: the displacementLaplacian solve
//  (the structure's deformation and control-surface rotations, set on
//  pointDisplacement in the vehicle's reference frame) followed by the
//  vehicle's rigid-body pose, which moves and rotates the whole domain:
//      x = centre + displacement + R (points0 + d - centre)
//  The pose comes from the coupling (the adapter's BodyPose readers, see
//  FSI/BodyFrame.H); `centre` is the rotation centre in the reference
//  configuration (the vehicle's CG). The Laplacian treats every direction
//  alike, so solving it on the rigidly moved mesh gives the same d.
//
//  The whole domain translates with the vehicle, but only the region within
//  `rotationInner` of the centre rotates rigidly with it; the rotation fades
//  (cos^2) to none at `rotationOuter`. Rotating the whole domain would sweep
//  the far-field boundaries through the air at omega * r (metres per second
//  per rad/s, metres away), forcing large flows through them.
//
//  Usage (constant/dynamicMeshDict):
//      mover
//      {
//          type        bodyDisplacementLaplacian;
//          libs        ("libpreciceAdapterFunctionObject.so");
//          centre      (0.48 0.25 0.018);
//          rotationInner 0.35;   // rigid rotation within this radius
//          rotationOuter 0.75;   // no rotation beyond this radius
//          diffusivity inverseDistance (airframe);
//      }
class bodyDisplacementLaplacian : public displacementLaplacian
{
    //- Rotation centre in the reference configuration
    const vector centre_;

    //- Radii of rigid rotation and of no rotation
    const scalar rotationInner_;
    const scalar rotationOuter_;

    //- Fraction of the rotation each point takes (1 within rotationInner_,
    //  0 beyond rotationOuter_), from its reference position
    scalarField rotationWeight_;

public:
    //- Runtime type information
    TypeName("bodyDisplacementLaplacian");

    //- Construct from polyMesh and dictionary
    bodyDisplacementLaplacian(const polyMesh&, const dictionary&);

    //- Construct from fvMesh and dictionary
    bodyDisplacementLaplacian(fvMesh&, const dictionary&);

    //- Destructor
    ~bodyDisplacementLaplacian() = default;

    //- New point positions: the deformation, then the rigid-body pose
    virtual tmp<pointField> newPoints();
};

}
}

#endif
