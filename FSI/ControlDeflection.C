#include "ControlDeflection.H"
#include "Utilities.H"

#include "pointPatchField.H"
#include "fixedValuePointPatchField.H"

#include <algorithm>

using namespace Foam;

// Rotate a point p about the line through hinge with direction axis by angle
// theta (radians). Uses the Rodrigues rotation formula.
static Foam::vector rotateAboutAxis(
    const Foam::vector& p,
    const Foam::vector& hinge,
    const Foam::vector& axis,
    const scalar theta)
{
    const Foam::vector d = p - hinge;
    const Foam::vector a = normalised(axis);
    const scalar c = Foam::cos(theta);
    const scalar s = Foam::sin(theta);

    return hinge + c*d + s*(a ^ d) + (1.0 - c)*(a & d)*a;
}

preciceAdapter::FSI::ControlDeflection::ControlDeflection(
    const Foam::fvMesh& mesh,
    std::vector<ControlSurfaceConfig>& controlSurfaces,
    bool servoAngle)
: mesh_(mesh),
  controlSurfaces_(&controlSurfaces),
  servoAngle_(servoAngle)
{
    dataType_ = scalar;

    // Readers are created while the interfaces are configured, before any
    // data is read, so the command reader sees this flag from the first step.
    if (servoAngle_)
    {
        for (auto& surface : *controlSurfaces_)
        {
            surface.servoDriven = true;
        }
    }

    // Read the name of the pointDisplacement field (if different)
    const dictionary& FSIdict =
        mesh_.lookupObject<IOdictionary>("preciceDict").subOrEmptyDict("FSI");
    namePointDisplacement_ = FSIdict.lookupOrDefault<word>(
        "namePointDisplacement", "pointDisplacement");
}

void preciceAdapter::FSI::ControlDeflection::initialize()
{
    storeReferencePoints();
}

void preciceAdapter::FSI::ControlDeflection::storeReferencePoints()
{
    referencePoints_.clear();

    for (const auto& surface : *controlSurfaces_)
    {
        const label patchID = mesh_.boundary().findIndex(surface.patch);

        if (patchID == -1)
        {
            FatalErrorInFunction
                << "Control surface patch '" << surface.patch
                << "' does not exist."
                << exit(FatalError);
        }

        // Initial (reference) positions of the patch points.
        const pointField localPoints =
            mesh_.boundary()[patchID].poly().localPoints();

        std::vector<double> refs;
        refs.reserve(localPoints.size() * 3);
        for (const point& p : localPoints)
        {
            refs.push_back(p.x());
            refs.push_back(p.y());
            refs.push_back(p.z());
        }
        referencePoints_.push_back(std::move(refs));
    }

    findHingeLinks();
}

void preciceAdapter::FSI::ControlDeflection::findHingeLinks()
{
    hingeLinks_.assign(controlSurfaces_->size(), {});

    std::vector<label> surfacePatches;
    for (const auto& surface : *controlSurfaces_)
    {
        surfacePatches.push_back(mesh_.boundary().findIndex(surface.patch));
    }

    for (std::size_t s = 0; s < controlSurfaces_->size(); ++s)
    {
        const ControlSurfaceConfig& surface = controlSurfaces_->at(s);
        const Foam::vector hinge(surface.hinge[0], surface.hinge[1], surface.hinge[2]);
        const Foam::vector axis = normalised(Foam::vector(
            surface.axis[0], surface.axis[1], surface.axis[2]));

        const polyPatch& patch = mesh_.boundary()[surfacePatches[s]].poly();
        const pointField& localPoints = patch.localPoints();
        const labelList& meshPoints = patch.meshPoints();

        // The surface's points that another (non-control-surface) patch also
        // has are where it is hinged to that patch; their position along the
        // axis is what the hinge displacement is interpolated over.
        forAll(localPoints, i)
        {
            const Foam::vector d = localPoints[i] - hinge;
            forAll(mesh_.boundary(), q)
            {
                if (std::find(surfacePatches.begin(), surfacePatches.end(), q)
                    != surfacePatches.end())
                {
                    continue;
                }
                const auto& pointMap = mesh_.boundary()[q].poly().meshPointMap();
                const auto iter = pointMap.find(meshPoints[i]);
                if (iter != pointMap.end())
                {
                    hingeLinks_[s].push_back({d & axis, q, iter()});
                    break;
                }
            }
        }

        std::sort(hingeLinks_[s].begin(), hingeLinks_[s].end(),
                  [](const HingeLink& a, const HingeLink& b) {
                      return a.span < b.span;
                  });

        if (!hingeLinks_[s].empty())
        {
            adapterInfo(
                "Control surface " + surface.patch + " follows "
                    + mesh_.boundary()[hingeLinks_[s].front().patchID].name()
                    + " at " + std::to_string(hingeLinks_[s].size())
                    + " shared hinge point(s).",
                "info");
        }
    }
}

std::size_t preciceAdapter::FSI::ControlDeflection::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    notImplemented("Writing control-surface deflections is not implemented!");
    return 0;
}

void preciceAdapter::FSI::ControlDeflection::read(double* buffer, const unsigned int dim)
{
    // Look up the pointDisplacement field (the dynamic mesh solver consumes it).
    if (!mesh_.foundObject<pointVectorField>(namePointDisplacement_))
        return;

    pointVectorField& pointDisplacement = const_cast<pointVectorField&>(
        mesh_.lookupObject<pointVectorField>(namePointDisplacement_));

    // One deflection angle per control surface, in buffer order.
    for (std::size_t s = 0; s < controlSurfaces_->size(); ++s)
    {
        ControlSurfaceConfig& surface = controlSurfaces_->at(s);
        const label patchID = mesh_.boundary().findIndex(surface.patch);
        if (patchID == -1)
            continue;

        const Foam::scalar theta = buffer[s];
        if (!servoAngle_)
        {
            surface.lastDeflection = theta;

            // The servo angle moves this surface; only relay the command.
            if (surface.servoDriven)
            {
                continue;
            }
        }

        const Foam::vector hinge(surface.hinge[0], surface.hinge[1], surface.hinge[2]);
        const Foam::vector axis(surface.axis[0], surface.axis[1], surface.axis[2]);

        // Point boundary field of the control-surface patch.
        vectorField& pField =
            refCast<vectorField>(pointDisplacement.boundaryFieldRef()[patchID]);

        const std::vector<double>& refs = referencePoints_.at(s);

        // Displacement of the hinge line: that of the patch the surface is
        // hinged to at the shared hinge points (zero if there are none, or if
        // that patch is not moved by a prescribed displacement).
        const std::vector<HingeLink>& links = hingeLinks_.at(s);
        std::vector<double> span;
        std::vector<Foam::vector> hingeDisp;
        for (const HingeLink& link : links)
        {
            const auto& other = pointDisplacement.boundaryField()[link.patchID];
            if (isA<fixedValuePointPatchVectorField>(other))
            {
                span.push_back(link.span);
                hingeDisp.push_back(refCast<const vectorField>(other)[link.index]);
            }
        }
        const Foam::vector a = normalised(axis);
        auto hingeDisplacement = [&](const Foam::scalar sp) -> Foam::vector {
            if (span.empty())
            {
                return Foam::vector::zero;
            }
            if (sp <= span.front())
            {
                return hingeDisp.front();
            }
            if (sp >= span.back())
            {
                return hingeDisp.back();
            }
            const std::size_t k =
                std::upper_bound(span.begin(), span.end(), sp) - span.begin();
            const Foam::scalar w = (sp - span[k - 1]) / (span[k] - span[k - 1]);
            return (1 - w) * hingeDisp[k - 1] + w * hingeDisp[k];
        };

        // Displacement = rotated - reference, plus the hinge-line
        // displacement: the patch rigidly rotates about the hinge and moves
        // with it, so it stays attached to a deforming patch it is hinged to.
        forAll(pField, i)
        {
            const Foam::vector p0(refs[3*i + 0], refs[3*i + 1], refs[3*i + 2]);
            const Foam::vector pRot = rotateAboutAxis(p0, hinge, axis, theta);
            pField[i] = pRot - p0 + hingeDisplacement((p0 - hinge) & a);
        }
    }
}

bool preciceAdapter::FSI::ControlDeflection::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::fixedPoints);
}

std::string preciceAdapter::FSI::ControlDeflection::getDataName() const
{
    return servoAngle_ ? "ServoAngle" : "Deflection";
}