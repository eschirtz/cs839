
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
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;

    std::array<int, 2> m_cellSize; // dimensions in grid cells
    T m_gridDX;
    int m_nFrames;

    static constexpr int m_pinchRadius = 5;

    void initialize()
    {
        initializeUSD("demo.usda");

        // Create a Cartesian lattice topology
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
    }

    void prepareFrame(const int frame)
    {
      // Grab the side and wave it
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
        // and on top
        currID = gridToParticleID(node_i , m_cellSize[0] );
        m_particleX[currID].Set(
          m_particleX[currID].data()[0],
          m_particleX[currID].data()[1],
          -zHeight
        );
      }
    }

    void relaxFreeNodes(const int iterations)
    {
        // Relax every interior node
        for(int iteration = 0; iteration < iterations; iteration++)
            for(int node_i = 1; node_i < m_cellSize[0]; node_i++)
            for(int node_j = 1; node_j < m_cellSize[1]; node_j++){

                if( std::abs(node_i - m_cellSize[0]/2) <= m_pinchRadius &&
                    std::abs(node_j - m_cellSize[1]/2) <= m_pinchRadius ) continue;

                int pCenter = gridToParticleID(node_i  ,node_j  );
                int pPlusX  = gridToParticleID(node_i+1,node_j  );
                int pMinusX = gridToParticleID(node_i-1,node_j  );
                int pPlusY  = gridToParticleID(node_i  ,node_j+1);
                int pMinusY = gridToParticleID(node_i  ,node_j-1);

                // Set the particle to be the average of it's neighbors
                m_particleX[pCenter] = .25 * ( m_particleX[pPlusX] + m_particleX[pMinusX] + m_particleX[pPlusY] + m_particleX[pMinusY]);
            }
    }

    void simulateFrame(const int frame)
    {
        // Relax every interior node
        for(int iteration = 0; iteration < 1; iteration++)
            for(int node_i = 1; node_i < m_cellSize[0]; node_i++)
            for(int node_j = 1; node_j < m_cellSize[1]; node_j++){

                int pCenter = gridToParticleID(node_i  ,node_j  );
                int pPlusX  = gridToParticleID(node_i+1,node_j  );
                int pMinusX = gridToParticleID(node_i-1,node_j  );
                int pPlusY  = gridToParticleID(node_i  ,node_j+1);
                int pMinusY = gridToParticleID(node_i  ,node_j-1);

                m_particleX[pCenter] = .25 * ( m_particleX[pPlusX] + m_particleX[pMinusX] + m_particleX[pPlusY] + m_particleX[pMinusY]);
            }
    }

private:
    inline int gridToParticleID(const int i, const int j) { return i * (m_cellSize[1]+1) + j; }
};

int main(int argc, char *argv[])
{
    LatticeMesh<float> simulationMesh;
    simulationMesh.m_cellSize = { 40, 40 };
    simulationMesh.m_gridDX = 0.025;
    simulationMesh.m_nFrames = 400;

    // Initialize the simulation example
    simulationMesh.initialize();
    simulationMesh.relaxFreeNodes(1000);

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
