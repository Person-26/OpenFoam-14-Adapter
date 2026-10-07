#include "BodyPose.H"
#include "BodyFrame.H"

using namespace Foam;

preciceAdapter::FSI::BodyPose::BodyPose(const Foam::fvMesh& mesh, bool rotation)
: mesh_(mesh),
  fieldName_(rotation ? "bodyRotation" : "bodyDisplacement")
{
    dataType_ = vector;
    bodyPoseField(mesh_, fieldName_);

    if (rotation)
    {
        const dictionary& FSIdict =
            mesh_.lookupObject<IOdictionary>("preciceDict").subOrEmptyDict("FSI");
        wallPatches_ =
            FSIdict.lookupOrDefault<wordList>("bodyWallPatches", wordList());
    }
}

void preciceAdapter::FSI::BodyPose::setWallVelocity() const
{
    // The walls move with the vehicle's rigid-body velocity, v + omega x r,
    // from the smooth pose curve (see read()). The plate's bending moves the
    // mesh but is deliberately left out of the wall velocity: a foam plate
    // this light has an added mass comparable to its own, and feeding its
    // bending velocity back through the explicit Fluid-Solid coupling (as
    // movingWallVelocity would) is unstable.
    if (wallPatches_.empty())
    {
        return;
    }
    const auto rate = [&](const word& name)
    {
        const Foam::scalar tau =
            mesh_.time().value() - bodyPoseTime(mesh_, name).value();
        return bodyPoseField(mesh_, name + "Rate").value()
             + tau * (2 * bodyPoseField(mesh_, name + "Acc").value()
             + tau * 3 * bodyPoseField(mesh_, name + "Jerk").value());
    };
    const Foam::vector v = rate("bodyDisplacement");
    // Rotation-vector rate: the angular velocity to first order in the angle
    // (the vehicle's attitude changes stay moderate here).
    const Foam::vector omega = rate("bodyRotation");
    const RigidPose pose = bodyPose(mesh_);
    const Foam::vector centre = pose.centre + pose.displacement;

    volVectorField& U =
        const_cast<volVectorField&>(mesh_.lookupObject<volVectorField>("U"));
    for (const word& patchName : wallPatches_)
    {
        const label patchID = mesh_.boundary().findIndex(patchName);
        if (patchID == -1)
        {
            continue;
        }
        const vectorField& Cf = mesh_.boundary()[patchID].Cf();
        vectorField& Ub = U.boundaryFieldRef()[patchID];
        forAll(Ub, i)
        {
            Ub[i] = v + (omega ^ (Cf[i] - centre));
        }
    }
}

std::size_t preciceAdapter::FSI::BodyPose::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    notImplemented("Writing the body pose is not implemented!");
    return 0;
}

void preciceAdapter::FSI::BodyPose::read(double* buffer, const unsigned int dim)
{
    // One vector on the single fixed point (every rank has it, see
    // Interface::readCouplingData()).
    Foam::vector v(Zero);
    for (unsigned int d = 0; d < dim; ++d)
    {
        v[d] = buffer[d];
    }
    // A new window's sample differs from the latest (the readers are called
    // every solver step, with the same value within a window). Keep the last
    // four samples. Over the window after a sample arrives, the pose follows
    // the Catmull-Rom segment between the second- and third-latest ones:
    // continuous in position and velocity, so the walls never jump in speed
    // (an added-mass pressure spike each window), two windows behind the
    // vehicle. bodyPose() evaluates the cubic in the time since arrival.
    const Foam::scalar t = mesh_.time().value();
    if (!haveSample_)
    {
        for (auto& p : samples_)
        {
            p = v;
        }
        tLatest_ = t;
        window_ = 0;
        haveSample_ = true;
    }
    else if (mag(v - samples_[3]) > 1e-15)
    {
        window_ = t - tLatest_;
        tLatest_ = t;
        for (int k = 0; k < 3; ++k)
        {
            samples_[k] = samples_[k + 1];
        }
        samples_[3] = v;
    }

    // Hermite segment from P1 = samples_[1] to P2 = samples_[2] over one
    // window W, with central-difference tangents, as a cubic in tau.
    const Foam::vector& P0 = samples_[0];
    const Foam::vector& P1 = samples_[1];
    const Foam::vector& P2 = samples_[2];
    const Foam::vector& P3 = samples_[3];
    Foam::vector a1(Zero), a2(Zero), a3(Zero);
    if (window_ > 0)
    {
        const Foam::scalar W = window_;
        const Foam::vector m1 = (P2 - P0) / (2 * W);
        const Foam::vector m2 = (P3 - P1) / (2 * W);
        a1 = m1;
        a2 = (3 * (P2 - P1) / W - 2 * m1 - m2) / W;
        a3 = (2 * (P1 - P2) / W + m1 + m2) / (W * W);
    }

    bodyPoseField(mesh_, fieldName_).value() = P1;
    bodyPoseField(mesh_, fieldName_ + "Rate").value() = a1;
    bodyPoseField(mesh_, fieldName_ + "Acc").value() = a2;
    bodyPoseField(mesh_, fieldName_ + "Jerk").value() = a3;
    bodyPoseTime(mesh_, fieldName_).value() = tLatest_;

    setWallVelocity();
}

bool preciceAdapter::FSI::BodyPose::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::fixedPoints);
}

std::string preciceAdapter::FSI::BodyPose::getDataName() const
{
    return fieldName_ == "bodyRotation" ? "BodyRotation" : "BodyDisplacement";
}
