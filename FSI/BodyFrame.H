#ifndef FSI_BODYFRAME_H
#define FSI_BODYFRAME_H

#include "objectRegistry.H"
#include "Time.H"
#include "tensor.H"
#include "uniformDimensionedFields.H"

namespace preciceAdapter
{
namespace FSI
{

//- Pose of a vehicle whose fluid mesh moves rigidly with it (the
//  bodyDisplacementLaplacian mesh mover). It is kept in the mesh's registry
//  as three uniformDimensionedVectorFields: "bodyCentre" (the rotation
//  centre, in the reference configuration), "bodyDisplacement" (the
//  translation of that centre) and "bodyRotation" (the rotation vector,
//  axis * angle). A reference point p is at
//      centre + displacement + R (p - centre)
//  now. Without the fields the pose is the identity (a fixed mesh).
//
//  The coupling delivers the pose once per window, at the window's start; a
//  mesh moved only then would jump every window (a spike of wall velocity,
//  with the added-mass pressure that goes with it). So the readers store a
//  cubic in the time since the latest sample arrived ("bodyPoseTime"): the
//  value and "<field>Rate", "<field>Acc", "<field>Jerk" coefficients of a
//  Catmull-Rom segment through the samples (see BodyPose::read()),
//  continuous in position and velocity, two windows behind the vehicle.
struct RigidPose
{
    Foam::vector centre = Foam::Zero;
    Foam::vector displacement = Foam::Zero;
    Foam::tensor R = Foam::tensor::I;

    //- Current position of a reference point
    Foam::vector toCurrent(const Foam::vector& p) const
    {
        return centre + displacement + (R & (p - centre));
    }

    //- Reference position of a current point
    Foam::vector toReference(const Foam::vector& p) const
    {
        return centre + (R.T() & (p - centre - displacement));
    }

    //- A current (inertial) vector in the body (reference) frame
    Foam::vector toBody(const Foam::vector& v) const
    {
        return R.T() & v;
    }
};

//- Rotation tensor of a rotation vector (Rodrigues)
inline Foam::tensor rotationTensor(const Foam::vector& w)
{
    const Foam::scalar theta = Foam::mag(w);
    if (theta < Foam::small)
    {
        return Foam::tensor::I;
    }
    const Foam::vector k = w / theta;
    const Foam::tensor K(
        0, -k.z(), k.y(),
        k.z(), 0, -k.x(),
        -k.y(), k.x(), 0);
    return Foam::tensor::I * Foam::cos(theta) + K * Foam::sin(theta)
         + (1 - Foam::cos(theta)) * (k * k);
}

//- Rotation vector of a rotation tensor (inverse of rotationTensor)
inline Foam::vector rotationVector(const Foam::tensor& R)
{
    const Foam::scalar c =
        Foam::max(-1.0, Foam::min(1.0, 0.5 * (Foam::tr(R) - 1)));
    const Foam::scalar theta = Foam::acos(c);
    const Foam::vector v(R.zy() - R.yz(), R.xz() - R.zx(), R.yx() - R.xy());
    const Foam::scalar s = Foam::sin(theta);
    if (s > 1e-6)
    {
        return v * (theta / (2 * s));
    }
    return 0.5 * v; // small angle (the vehicle never turns near 180 degrees here)
}

//- The named pose field, created (with value `init`) if it does not exist
inline Foam::uniformDimensionedVectorField& bodyPoseField(
    const Foam::objectRegistry& db,
    const Foam::word& name,
    const Foam::vector& init = Foam::Zero)
{
    if (!db.foundObject<Foam::uniformDimensionedVectorField>(name))
    {
        Foam::uniformDimensionedVectorField* f =
            new Foam::uniformDimensionedVectorField(
                Foam::IOobject(
                    name,
                    db.time().name(),
                    "uniform",
                    db,
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::NO_WRITE),
                Foam::dimensionedVector(name, Foam::dimless, init));
        f->store();
    }
    return const_cast<Foam::uniformDimensionedVectorField&>(
        db.lookupObject<Foam::uniformDimensionedVectorField>(name));
}

//- The pose sample time field, created (zero) if it does not exist
inline Foam::uniformDimensionedScalarField& bodyPoseTime(
    const Foam::objectRegistry& db)
{
    const Foam::word name("bodyPoseTime");
    if (!db.foundObject<Foam::uniformDimensionedScalarField>(name))
    {
        Foam::uniformDimensionedScalarField* f =
            new Foam::uniformDimensionedScalarField(
                Foam::IOobject(
                    name,
                    db.time().name(),
                    "uniform",
                    db,
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::NO_WRITE),
                Foam::dimensionedScalar(name, Foam::dimTime, 0));
        f->store();
    }
    return const_cast<Foam::uniformDimensionedScalarField&>(
        db.lookupObject<Foam::uniformDimensionedScalarField>(name));
}

//- The pose at the current time (the identity if there is none)
inline RigidPose bodyPose(const Foam::objectRegistry& db)
{
    RigidPose pose;
    if (db.foundObject<Foam::uniformDimensionedVectorField>("bodyCentre"))
    {
        const Foam::scalar tau =
            db.time().value() - bodyPoseTime(db).value();
        pose.centre = bodyPoseField(db, "bodyCentre").value();
        const auto cubic = [&](const Foam::word& name)
        {
            return bodyPoseField(db, name).value()
                 + tau * (bodyPoseField(db, name + "Rate").value()
                 + tau * (bodyPoseField(db, name + "Acc").value()
                 + tau * bodyPoseField(db, name + "Jerk").value()));
        };
        pose.displacement = cubic("bodyDisplacement");
        pose.R = rotationTensor(cubic("bodyRotation"));
    }
    return pose;
}

}
}

#endif
