
#include "pxr/pxr.h"


#include "pxr/usd/sdf/layer.h"
#include "pxr/usd/sdf/path.h"
#include "pxr/usd/usd/stage.h"
#include "pxr/usd/usdGeom/mesh.h"
#include "pxr/base/vt/array.h"

#include "pxr/base/gf/range3f.h"

#include <iostream>
#include <cmath>

PXR_NAMESPACE_USING_DIRECTIVE

template<class T, int d> // d is the dimension of the mesh elements, e.g. 3 for triangles, 4 for quads
struct AnimatedMesh
{
    SdfLayerRefPtr m_layer;
    UsdStageRefPtr m_stage;
    UsdGeomMesh m_mesh;
    UsdAttribute m_pointsAttribute;

    GfRange3f m_extent;
    int m_lastFrame;

    std::vector<std::array<int, d>> m_meshElements;
    std::vector<GfVec3f> m_particleX;

    AnimatedMesh()
        :m_lastFrame(-1)
    {}

    void initializeUSD(const std::string filename)
    {
        // Create the layer to populate.
        m_layer = SdfLayer::CreateNew(filename);

        // Create a UsdStage with that root layer.
        m_stage = UsdStage::Open(m_layer);
    }

    void initializeTopology()
    {
        // Create a mesh for this surface
        m_mesh = UsdGeomMesh::Define(m_stage, SdfPath("/MeshSurface"));

        // Create appropriate buffers for vertex counts and indices, and populate them
        VtIntArray faceVertexCounts, faceVertexIndices;
        for (const auto& element : m_meshElements) {
            faceVertexCounts.push_back(element.size());
            for (const auto& vertex : element)
                faceVertexIndices.push_back(vertex);
        }

        // Now set the attributes
        m_mesh.GetFaceVertexCountsAttr().Set(faceVertexCounts);
        m_mesh.GetFaceVertexIndicesAttr().Set(faceVertexIndices);
    }

    void initializeParticles()
    {
        // Grab the points (Positions) attribute, and indicate it is time-varying
        m_pointsAttribute = m_mesh.GetPointsAttr();
        m_pointsAttribute.SetVariability(SdfVariabilityVarying);
    }

    void writeFrame(const int frame)
    {
        std::cout << "Writing frame " << frame << " ..." << std::endl;

        // Check that there are any particles to write at all
        if (m_particleX.empty())
            throw std::logic_error("Empty array of input vertices");

        // Check that frames have been written in sequence
        if(frame != m_lastFrame+1)
            throw std::logic_error("Non-consequtive frame sequence requested in writeFrame()");
        m_lastFrame = frame;

        // Update extent
        for (const auto& pt : m_particleX)
            m_extent.UnionWith(pt);

        // Copy particleX into VtVec3fArray for Usd
        VtVec3fArray usdPoints;
        usdPoints.assign(m_particleX.begin(), m_particleX.end());

        // Write the points attribute for the given frame
        m_pointsAttribute.Set(usdPoints, (double) frame);
    }

    void writeUSD()
    {
        // Set up the timecode
        m_stage->SetStartTimeCode(0.);
        m_stage->SetEndTimeCode((double) m_lastFrame);

        // Set the effective extent
        VtVec3fArray extentArray(2);
        extentArray[0] = m_extent.GetMin();
        extentArray[1] = m_extent.GetMax();
        m_mesh.GetExtentAttr().Set(extentArray);

        // Save USD file
        m_stage->GetRootLayer()->Save();
        std::cout << "USD file saved!" << std::endl;
    }
};

template<class T>
struct LatticeMesh : public AnimatedMesh<T, 4>
{
    using Base = AnimatedMesh<T, 4>;
    using Base::m_meshElements;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Base::m_particleX; // positional data array

    std::vector<GfVec3f> m_particleV; // velocity data array
    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;
    // Simulation
    int m_nFrames;
    int m_subSteps;
    T m_timeValuePerFrame = 1;
    // Physics
    T m_particleMass = 0.90;
    T m_stiffnessCoeff = 0.75;
    T m_dampingCoeff = 0.2;

    static constexpr int m_pinchRadius = 5;

    void initialize()
    {
        initializeUSD("demo.usda");

        // Create a Cartesian lattice topology ( add the elements )
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++)
        for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++)
            m_meshElements.emplace_back(
                std::array<int, 4>{
                    gridToParticleID(cell_i  , cell_j  ),
                    gridToParticleID(cell_i+1, cell_j  ),
                    gridToParticleID(cell_i+1, cell_j+1),
                    gridToParticleID(cell_i  , cell_j+1)
                }
            );
        initializeTopology();

        // Also initialize the associated particles
        for(int node_i = 0; node_i <= m_cellSize[0]; node_i++)
        for(int node_j = 0; node_j <= m_cellSize[1]; node_j++)
            // pushing points into the array
            m_particleX.emplace_back(m_gridDX * (T)node_i, m_gridDX * (T)node_j, T());
        initializeParticles();

        // Also resize the velocity array to match
        m_particleV.resize(m_particleX.size());
    }

    void prepareFrame(const int frame)
    {
      // Grab the side and wave it sinusoidally
      float height = 0.5;
      float speed = 0.125;
      float zHeight = height * sin((T)frame * speed);
      for(int node_i = 0; node_i <= m_cellSize[0]; node_i++){
        int currID = gridToParticleID(node_i  , 0 );
        // Update all on bottom
        m_particleX[currID].Set(
          m_particleX[currID].data()[0],
          m_particleX[currID].data()[1],
          zHeight
        );
      }
    }

    void setForceValues(std::vector<GfVec3f>& f)
    {
      // Must initialize every single point! Remove once fix large for loop
      for(int i=0; i<f.size(); i++)
        f[i] = GfVec3f(0.0,0.0,0.0);

      for(int node_i = 0; node_i <= m_cellSize[0]; node_i++)
      for(int node_j = 0; node_j <= m_cellSize[1]; node_j++){
        // Get indices of neighbors in the mesh
        int pCenter = gridToParticleID(node_i  ,node_j  );
        int pPlusX  = node_i+1 > m_cellSize[0] ? pCenter : gridToParticleID(node_i+1,node_j);
        int pMinusX = node_i-1 < 0 ? pCenter : gridToParticleID(node_i-1,node_j);
        int pPlusY  = node_j+1 > m_cellSize[1] ? pCenter : gridToParticleID(node_i  ,node_j+1);
        int pMinusY = node_j-1 < 0 ? pCenter : gridToParticleID(node_i  ,node_j-1);

        GfVec3f diffVec;
        GfVec3f restVec;
        // Positive-X neighbor spring force
        diffVec = (m_particleX[pCenter] - m_particleX[pPlusX]);
        restVec = diffVec.GetNormalized() * m_gridDX;
        f[pCenter] -= m_stiffnessCoeff * (diffVec - restVec);
        f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pPlusX]);
        // Negative-X neighbor spring force
        diffVec = (m_particleX[pCenter] - m_particleX[pMinusX]);
        restVec = diffVec.GetNormalized() * m_gridDX;
        f[pCenter] -= m_stiffnessCoeff * (diffVec - restVec);
        f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pMinusX]);
        // Positive-Y neighbor spring force
        diffVec = (m_particleX[pCenter] - m_particleX[pPlusY]);
        restVec = diffVec.GetNormalized() * m_gridDX;
        f[pCenter] -= m_stiffnessCoeff * (diffVec - restVec);
        f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pPlusY]);
        // Negative-Y neighbor spring force
        diffVec = (m_particleX[pCenter] - m_particleX[pMinusY]);
        restVec = diffVec.GetNormalized() * m_gridDX;
        f[pCenter] -= m_stiffnessCoeff * (diffVec - restVec);
        f[pCenter] -= m_dampingCoeff * (m_particleV[pCenter] - m_particleV[pMinusY]);
      }
    }

    void simulateSubstep(const T dt)
    {
        const int nParticles = m_particleX.size();
        std::vector<GfVec3f> forceVector(nParticles);
        setForceValues(forceVector);

        for(int p = 0; p < nParticles; p++){
          // a=F/m, a * t = v. Add the change in velocity
          m_particleV[p] = (dt / m_particleMass ) * forceVector[p]; //(dt / m_particleMass) * forceVector[p];
        }
        for(int p = 0; p < nParticles; p++){
          // Update position
          m_particleX[p] += dt * m_particleV[p];
        }
    }

    void simulateFrame(const int frame)
    {
      T dT = m_timeValuePerFrame / (T) m_subSteps;
      for(int step = 0; step < m_subSteps; step++)
        simulateSubstep(dT);
    }

private:
    inline int gridToParticleID(const int i, const int j) { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 30, 30 };
    simulationMesh.m_gridDX = 0.025;
    simulationMesh.m_nFrames = 400;
    simulationMesh.m_subSteps = 5;

    // Initialize the simulation example
    simulationMesh.initialize();

    // Output the initial shape of the mesh
    simulationMesh.writeFrame(0);

    // Perform the animation, output results at each frame
    for(int frame = 1; frame <= simulationMesh.m_nFrames; frame++){
        simulationMesh.prepareFrame(frame);
        simulationMesh.simulateFrame(frame);
        simulationMesh.writeFrame(frame);
    }

    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}
