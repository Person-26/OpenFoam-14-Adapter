#include "PropellerLoad.H"
#include "BodyFrame.H"

#include "fvModels.H"
#include "propellerDisk.H"
#include "rotorDisk.H"
#include "vector.H"

using namespace Foam;

preciceAdapter::FSI::PropellerLoad::PropellerLoad(
    const Foam::fvMesh& mesh,
    std::vector<std::string> propellerNames,
    bool isMoment,
    Foam::scalar rho)
: mesh_(mesh),
  propellerNames_(std::move(propellerNames)),
  isMoment_(isMoment),
  rho_(rho)
{
    dataType_ = vector;
}

std::size_t preciceAdapter::FSI::PropellerLoad::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    // Get the fvModels object (contains all the fvModels of the case)
    const fvModels& models = mesh_.lookupObject<fvModels>("fvModels");

    std::size_t bufferIndex = 0;

    // Write one value (thrust or torque) per propeller hub, in the same
    // order as the propeller fvModels were specified.
    for (const std::string& name : propellerNames_)
    {
        // Either model gives the fluid's load on the propeller (the thrust
        // and reaction torque on the airframe as they are; propellerDisk's
        // normal is the thrust direction, it pushes the fluid along -normal),
        // times rho for an incompressible solver's kinematic loads. A
        // blade-element rotorDisk also gives the in-plane force and the hub
        // moments of oblique inflow. In the body frame when the mesh moves
        // with a vehicle.
        const fvModel& model = models[name];
        Foam::vector load;
        if (isA<fv::rotorDisk>(model))
        {
            const fv::rotorDisk& rotor = refCast<const fv::rotorDisk>(model);
            load = isMoment_ ? rotor.moment() : rotor.force();
        }
        else
        {
            const fv::propellerDisk& prop = dynamicCast<const fv::propellerDisk>(model);
            load = isMoment_ ? prop.moment() : prop.force();
        }
        const Foam::vector value = bodyPose(mesh_).toBody(rho_ * load);

        for (unsigned int d = 0; d < dim; ++d)
            buffer[bufferIndex++] = value[d];
    }

    return bufferIndex;
}

void preciceAdapter::FSI::PropellerLoad::read(double* buffer, const unsigned int dim)
{
    notImplemented("Reading propeller loads is not implemented!");
}

bool preciceAdapter::FSI::PropellerLoad::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::fixedPoints);
}

std::string preciceAdapter::FSI::PropellerLoad::getDataName() const
{
    return isMoment_ ? "PropTorque" : "Thrust";
}