#include "PropellerRPM.H"

#include "fvModels.H"
#include "propellerDisk.H"

using namespace Foam;

preciceAdapter::FSI::PropellerRPM::PropellerRPM(
    const Foam::fvMesh& mesh,
    std::vector<std::string> propellerNames)
: mesh_(mesh),
  propellerNames_(std::move(propellerNames))
{
    dataType_ = scalar;
}

std::size_t preciceAdapter::FSI::PropellerRPM::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    notImplemented("Writing propeller RPM is not implemented!");
    return 0;
}

void preciceAdapter::FSI::PropellerRPM::read(double* buffer, const unsigned int dim)
{
    // Get the fvModels object (contains all the fvModels of the case).
    fvModels& models = mesh_.lookupObjectRef<fvModels>("fvModels");

    // One commanded RPM per propeller hub, in the same order as the
    // propeller fvModels were specified. The sign selects the rotation
    // direction (handled by propellerDisk::setRotationSpeed).
    for (std::size_t i = 0; i < propellerNames_.size(); ++i)
    {
        fv::propellerDisk& prop =
            dynamicCast<fv::propellerDisk>(models[propellerNames_[i]]);

        prop.setRotationSpeed(buffer[i]);
    }
}

bool preciceAdapter::FSI::PropellerRPM::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::fixedPoints);
}

std::string preciceAdapter::FSI::PropellerRPM::getDataName() const
{
    return "RPM";
}
