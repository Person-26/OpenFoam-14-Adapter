#include "ControlCommand.H"

using namespace Foam;

preciceAdapter::FSI::ControlCommand::ControlCommand(
    const Foam::fvMesh& mesh,
    std::vector<ControlSurfaceConfig>& controlSurfaces)
: mesh_(mesh),
  controlSurfaces_(&controlSurfaces)
{
    dataType_ = scalar;
}

void preciceAdapter::FSI::ControlCommand::initialize()
{
    // Nothing to do: the command is taken from the shared configuration.
}

std::size_t preciceAdapter::FSI::ControlCommand::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    // Relay the deflection that the ControlDeflection reader last applied: the
    // same command the flight controller sent, for the structural servo to
    // track.
    for (std::size_t s = 0; s < controlSurfaces_->size(); ++s)
        buffer[s] = controlSurfaces_->at(s).lastDeflection;

    return controlSurfaces_->size();
}

void preciceAdapter::FSI::ControlCommand::read(double* buffer, const unsigned int dim)
{
    notImplemented("Reading control-surface commands is not implemented!");
}

bool preciceAdapter::FSI::ControlCommand::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (meshConnectivity == false);
}

std::string preciceAdapter::FSI::ControlCommand::getDataName() const
{
    return "ServoCommand";
}
