#include "BodyDisplacementLaplacian.H"
#include "BodyFrame.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fvMotionSolvers
{
    defineTypeNameAndDebug(bodyDisplacementLaplacian, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        bodyDisplacementLaplacian,
        fvMesh
    );

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        bodyDisplacementLaplacian,
        dictionary
    );
}
}


Foam::fvMotionSolvers::bodyDisplacementLaplacian::bodyDisplacementLaplacian
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementLaplacian(mesh, dict),
    centre_(dict.lookup<vector>("centre")),
    rotationInner_(dict.lookupOrDefault<scalar>("rotationInner", great)),
    rotationOuter_(dict.lookupOrDefault<scalar>("rotationOuter", great)),
    rotationWeight_(points0().size(), 1)
{
    forAll(rotationWeight_, i)
    {
        const scalar r = mag(points0()[i] - centre_);
        if (r >= rotationOuter_)
        {
            rotationWeight_[i] = 0;
        }
        else if (r > rotationInner_)
        {
            const scalar s =
                (r - rotationInner_)/(rotationOuter_ - rotationInner_);
            rotationWeight_[i] = sqr(cos(0.5*constant::mathematical::pi*s));
        }
    }

    // Create the pose (the identity until the coupling sets it). The centre
    // lives here, so that the adapter's load writers and the propeller disks
    // use the same one.
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyCentre", centre_).value() =
        centre_;
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyDisplacement");
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyRotation");
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyDisplacementRate");
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyRotationRate");
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyDisplacementAcc");
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyRotationAcc");
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyDisplacementJerk");
    preciceAdapter::FSI::bodyPoseField(mesh, "bodyRotationJerk");
    preciceAdapter::FSI::bodyPoseTime(mesh, "bodyDisplacement");
    preciceAdapter::FSI::bodyPoseTime(mesh, "bodyRotation");
}


Foam::fvMotionSolvers::bodyDisplacementLaplacian::bodyDisplacementLaplacian
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    bodyDisplacementLaplacian(mesh.poly(), dict)
{}


Foam::tmp<Foam::pointField>
Foam::fvMotionSolvers::bodyDisplacementLaplacian::newPoints()
{
    // Reference configuration plus deformation (points0 + d)
    const pointField deformed(displacementLaplacian::newPoints());

    const preciceAdapter::FSI::RigidPose pose =
        preciceAdapter::FSI::bodyPose(fvMotionSolver::mesh());
    const vector w = preciceAdapter::FSI::rotationVector(pose.R);

    tmp<pointField> tnewPoints(new pointField(deformed.size()));
    pointField& newPoints = tnewPoints.ref();
    forAll(deformed, i)
    {
        // Full translation; the rotation scaled by the point's weight (the
        // vehicle region rotates rigidly, the far field not at all).
        const scalar wi = rotationWeight_[i];
        const tensor Ri =
            wi >= 1 ? pose.R
          : wi <= 0 ? tensor::I
          : preciceAdapter::FSI::rotationTensor(wi*w);
        newPoints[i] =
            pose.centre + pose.displacement + (Ri & (deformed[i] - pose.centre));
    }
    return tnewPoints;
}
