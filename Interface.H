#ifndef INTERFACE_H
#define INTERFACE_H

#include <map>
#include <string>
#include <vector>
#include "OpenFOAMCompatibility.H"
#include "CouplingDataUser.H"
#include <precice/precice.hpp>

#include "pointPatchField.H"

namespace preciceAdapter
{

class Interface
{
protected:
    //- preCICE solver interface
    precice::Participant& precice_;

    //- Mesh name used in the preCICE configuration
    std::string meshName_;

    //- Point locations type (faceCenters vs faceNodes)
    LocationType locationType_ {LocationType::none};

    //- Mesh ID assigned by preCICE to the interface
    int meshID_;

    //- Names of the OpenFOAM patches that form the interface
    std::vector<std::string> patchNames_;

    //- OpenFOAM patches that form the interface
    std::vector<int> patchIDs_;

    //- Names of the OpenFOAM cell sets to be coupled (for volume coupling)
    std::vector<std::string> cellSetNames_;

    //- Names of the propellerDisk fvModels to be coupled (for fixedPoints)
    std::vector<std::string> propellerNames_;

    //- Explicit fixed 3D points of the interface (for the fixedPoints
    //  location type, e.g. control-surface hinges)
    std::vector<std::vector<double>> interfacePoints_;

    //- Fixed 3D points of the interface (for the fixedPoints location type)
    std::vector<std::vector<double>> fixedPoints_;

    //- Number of data points (cell centers) on the interface
    int numDataLocations_ = 0;

    //- Vertex IDs assigned by preCICE
    std::vector<int> vertexIDs_;

    //- Buffer for the coupling data
    std::vector<double> dataBuffer_;

    //- Vector of CouplingDataReaders
    std::vector<CouplingDataUser*> couplingDataReaders_;

    //- Vector of CouplingDataWriters
    std::vector<CouplingDataUser*> couplingDataWriters_;

    //Switch for faceTriangulation (nearest projection)
    bool meshConnectivity_;

    //- Reset the displacement during interface definition
    bool restartFromDeformed_;

    // Simulation dimension
    unsigned int dim_;

    //- History of a predicted data field (see predict())
    struct History
    {
        long window = -1;
        std::vector<double> cur, prev, vcur, vprev;
    };

    //- Predictor mode, its velocity data, and the current window
    //  (set by the Adapter before each read)
    int predictor_ = 0;
    bool restarted_ = false;
    std::string velocityData_;
    long window_ = 0;
    double tau_ = 0;
    double windowSize_ = 0;
    std::map<std::string, History> history_;

    //- Replace the first n buffer entries (data held constant across the
    //  window by a parallel coupling) by the predictor's value at tau_
    void predict(const std::string& dataName, std::size_t n, double relativeReadTime);

    //- Extracts the locations of the face centers or face nodes
    //  and exposes them to preCICE with setMeshVertices
    void configureMesh(const Foam::fvMesh& mesh,
                       const std::string& namePointDisplacement,
                       const std::string& nameCellDisplacement);

public:
    //- Constructor
    Interface(
        precice::Participant& precice,
        const Foam::fvMesh& mesh,
        std::string meshName,
        std::string locationsType,
        std::vector<std::string> patchNames,
        std::vector<std::string> cellSetNames,
        std::vector<std::string> propellerNames,
        std::vector<std::vector<double>> interfacePoints,
        bool meshConnectivity,
        bool restartFromDeformed,
        const std::string& namePointDisplacement,
        const std::string& nameCellDisplacement);

    //- Add a CouplingDataUser to read data from the interface
    void addCouplingDataReader(
        std::string dataName,
        CouplingDataUser* couplingDataReader);


    //- Add a CouplingDataUser to write data on the interface
    void addCouplingDataWriter(
        std::string dataName,
        CouplingDataUser* couplingDataWriter);

    //- Allocate an appropriate buffer for scalar or vector data.
    //  If at least one couplingDataUser has vector data, then
    //  define a buffer for 3D data. Otherwise, for 1D data.
    void createBuffer();

    //- Predictors for the Displacement and ServoAngle a parallel coupling
    //  holds constant across a window (the latest sample, from the end of
    //  the previous window): with d the latest sample, dp the one before,
    //  v and vp the matching PlateVelocity samples, h the window size and
    //  tau the fraction of the window elapsed at the read time,
    //    hold      d
    //    ramp      dp + tau (d - dp)       (continuous, one window late)
    //    linear    d + tau (d - dp)
    //    velocity  d + tau h v
    //    ab2       d + tau h v + tau^2 h (v - vp) / 2
    //  velocity and ab2 need velocityData; without it (and for ServoAngle,
    //  which has no rate) they fall back to linear.
    enum PredictorMode { hold = 0, ramp, linear, velocity, ab2 };
    void setPredictor(int mode, const std::string& velocityData);

    //- The fluid restarted from a saved time: the first window's ServoAngle
    //  and RPM are preCICE's uninitialised zeros (a received mesh cannot
    //  carry initial data), so keep the saved elevon deflection and the
    //  rotor speeds the restart set in fvModels for that window instead.
    void setRestarted(bool restarted) { restarted_ = restarted; }
    void setWindow(long window, double tau, double windowSize);

    //- Call read() on each registered couplingDataReader to read the coupling
    //  data from the buffer and apply the boundary conditions
    void readCouplingData(double relativeReadTime);

    //- Call write() on each registered couplingDataWriter to extract the boundary
    //  data and write them into the buffer
    void writeCouplingData();

    //- Destructor
    ~Interface();
};

}

#endif
