#include "Stress.H"

#include "Utilities.H"

#include <algorithm>
#include <string>

using namespace Foam;

preciceAdapter::FSI::Stress::Stress(
    const Foam::fvMesh& mesh,
    const std::string solverType)
: ForceBase(mesh, solverType)
{
    Stress_ = new volVectorField(
        IOobject(
            "Stress",
            mesh_.time().timeName(mesh_.time().value()),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector(
            "pdim",
            dimensionSet(1, -1, -2, 0, 0, 0, 0),
            Foam::vector::zero));
}

std::size_t preciceAdapter::FSI::Stress::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    const std::size_t n = this->writeToBuffer(buffer, *Stress_, dim);

    // A zero-thickness baffle puts both of its sides into the interface at the
    // same face centre, each carrying the absolute traction of its own side.
    // A structure behind the baffle is loaded by their sum (the pressure
    // difference across it), and a nearest-neighbour mapping could only pick
    // one side arbitrarily. Write the net traction to both sides instead.
    if (!bafflePairsBuilt_)
    {
        findBafflePairs();
    }
    for (const auto& pair : bafflePairs_)
    {
        for (unsigned int d = 0; d < dim; ++d)
        {
            const double net =
                buffer[pair.first * dim + d] + buffer[pair.second * dim + d];
            buffer[pair.first * dim + d] = net;
            buffer[pair.second * dim + d] = net;
        }
    }
    return n;
}

void preciceAdapter::FSI::Stress::findBafflePairs()
{
    // Pairing is topological, so it is found once, on the geometry at the
    // first write. Faces of the interface in buffer order.
    std::vector<Foam::vector> centres;
    std::vector<Foam::vector> normals;
    std::vector<double> sizes;
    for (const label patchID : patchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchID];
        const vectorField& Cf = patch.Cf();
        const tmp<vectorField> tnf = patch.nf();
        const scalarField& magSf = patch.magSf();
        forAll(Cf, i)
        {
            centres.push_back(Cf[i]);
            normals.push_back(tnf()[i]);
            sizes.push_back(Foam::sqrt(magSf[i]));
        }
    }

    // Sweep along x: candidates for a face lie within its tolerance in x.
    std::vector<std::size_t> order(centres.size());
    for (std::size_t i = 0; i < order.size(); ++i)
    {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        return centres[a].x() < centres[b].x();
    });

    std::vector<bool> paired(centres.size(), false);
    for (std::size_t k = 0; k < order.size(); ++k)
    {
        const std::size_t a = order[k];
        if (paired[a])
        {
            continue;
        }
        const double tol = 1e-6 * sizes[a];
        for (std::size_t m = k + 1; m < order.size(); ++m)
        {
            const std::size_t b = order[m];
            if (centres[b].x() - centres[a].x() > tol)
            {
                break;
            }
            if (!paired[b] && Foam::mag(centres[b] - centres[a]) <= tol
                && (normals[a] & normals[b]) < -0.5)
            {
                bafflePairs_.emplace_back(a, b);
                paired[a] = true;
                paired[b] = true;
                break;
            }
        }
    }
    bafflePairsBuilt_ = true;

    if (!bafflePairs_.empty())
    {
        adapterInfo(
            "Stress: writing the net traction on "
                + std::to_string(bafflePairs_.size())
                + " two-sided baffle face pair(s) of the interface.",
            "info");
    }
}

void preciceAdapter::FSI::Stress::read(double* buffer, const unsigned int dim)
{
    this->readFromBuffer(buffer);
}

bool preciceAdapter::FSI::Stress::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FSI::Stress::getDataName() const
{
    return "Stress";
}

Foam::tmp<Foam::vectorField> preciceAdapter::FSI::Stress::getFaceVectors(const unsigned int patchID) const
{
    // face normal vectors
    return mesh_.boundary()[patchID].nf();
}

preciceAdapter::FSI::Stress::~Stress()
{
    delete Stress_;
}
