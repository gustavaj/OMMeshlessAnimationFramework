#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/hash.hpp>

#include <unordered_map>

#include "LatticeUtility.h"
#include "Simulators.h"

/*
	In general, all local surface and patch orderings are left to right, top to bottom e.g.
		p00 - p10
	^ v	 |     | 
	|	p01 - p11
	 --> u
	Where the first number indicates u, and the second v.


	TODO list:
		Sphere:
		-Set all the top and bottom local surfaces to point to the same matrix.

		Simulators:
		-Fix RandomSphereSimulator
		-Add a RotationSimulator(maxAngle, axis)
		-Make it possible to use more than one simulator

		Pre-evaluation:
		-Evaluate local surfaces on the CPU. (pos and first partial derivatives)
		-In vulkan implementation, save results to buffer/texture
		-Use in TES to improve(?) evaluation time.
		-Look at memory usage with respect to large lattices.
		-Experiment with different buffer/texture dimensions.

		Pixel-accurate rendering:
		-Set tessellation factors in the TCS based on some error metric.
*/

namespace OML {

	using Vec2f = OpenMesh::Vec2f;
	using Vec3f = OpenMesh::Vec3f;
	using Col3 = OpenMesh::Vec3uc;

	/*
		Not in use, currently only one possible local surface, will probably not change.
	*/
	enum class LocalSurfaceType
	{
		Bezier3x3 = 0
	};

	/*
		The loci are located at the knot vectors(?). At every locus, there is one local surface.
		This is made up from the points:
		[m_controlPoints[controlPointIndex], m_controlPoints[controlPointIndex + controlPointCount]]
		The transformation matrix of the local surface is at m_matrices[matrixIndex]
		The normal is calculated as the average of the normal of all the surrounding faces.
		The normal is used for simulation by the NormalSinSimulator.
		The boundaryIndices maps a patch to an index into the m_boundaries vector.
		The boundaries contain info for how far the local surface should be evaluated for the given patch.
	*/
	struct Locus
	{
		Locus() : controlPointIndex(0), controlPointCount(0), matrixIndex(0) {}

		uint32_t controlPointIndex;
		uint32_t controlPointCount;
		uint32_t matrixIndex;
		OpenMesh::VertexHandle vh;
		Vec3f normal;
		std::unordered_map<OpenMesh::FaceHandle, uint32_t> boundaryIndices;
	};

	/*
		The patch contains four indices to the loci that makes up the patch.
	*/
	struct Patch
	{
		glm::vec3 color;
		std::array<size_t, 4> lociIndices;
		OpenMesh::FaceHandle fh;
	};

	/*
		Different B-Functions that can be used.
		B1Poly = "3x^2-2x^3"
		B2Polt = "6x^5-15x^4+10x^3"
		LERBS  = "1/(1+e^(1/x-1/(1-x)))"
	*/
	enum class BFunctionType {
		B1Poly = 0, B2Poly, LERBS
	};
	// The possible B-functions
	const std::vector<std::string> BFunctionNames = {
		"3x^2-2x^3", "6x^5-15x^4+10x^3", "1/(1+e^(1/x-1/(1-x)))"
	};

	// Custom traits passed to the OpenMesh::PolyMesh class.
	struct LatticeTraits : public OpenMesh::DefaultTraits
	{
		typedef Vec3f Point;

		// Vertices are colored by their valence, status is used for deleting vertices.
		VertexAttributes(OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);

		// Edges are colored depending on if they are on the boudnary or inside the mesh
		// Status used for deleting edges.
		EdgeAttributes(OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);

		// Status used for deleting halfedges
		HalfedgeAttributes(OpenMesh::Attributes::Status);

		// Status used for deleting faces.
		FaceAttributes(OpenMesh::Attributes::Status);
	};

	enum class LocusType {
		Unresolved = 0, Corner, Boundary, Inner, T, T_Terminal
	};

	// Different properties of the vertices
	namespace LatticeProperties {
		// Holds the value of the valence of a vertex. So it does not need to calculated more than once.
		static OpenMesh::VPropHandleT<size_t> VertexValence;
		// Holds the index into the m_loci vector
		static OpenMesh::VPropHandleT<uint32_t> LocusIndex;
		static OpenMesh::VPropHandleT<LocusType> Type;
	};

	// Colors used for the grid. The color of a vertex is indexed by its valence.
	// Vertices with a valence less than 2 are discarded in the fragment shader.
	const std::vector<Col3> LOCUS_VALENCE_COLOR = {
		Col3(255,255,255), Col3(255,255,255), Col3(255,0,0), Col3(0,255,0), Col3(0,0,255), Col3(0,255,255)
	};
	const Col3 BOUNDARY_EDGE_COLOR = Col3(0, 0, 0);
	const Col3 INNER_EDGE_COLOR = Col3(200, 200, 200);

	class Lattice : public OpenMesh::PolyMesh_ArrayKernelT<LatticeTraits>
	{
	public:
		// Framework
		Lattice();
		Lattice(std::string name);
		~Lattice();

		// Functions for creating patches
		/* Add patch in space */
		void addPatch(Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight);
		/* Add patch in xy-plane */
		void addPatch(Vec2f topLeft, Vec2f topRight, Vec2f bottomLeft, Vec2f bottomRight);
		/* Add patch in xy-plane */
		void addPatch(Vec2f topLeft, float width, float height);
		/* Add patch in space */
		void addPatch(Vec3f topLeft, float width, float height, float depth);
		/* Add grid in xy-plane */
		void addGrid(Vec2f topLeft, float width, float height, int rows, int cols);
		/* Add donut in xy-plane */
		void addDonut(Vec2f topLeft, float outerRadius, float innerRadius);
		/* Add cylinder in space */
		void addCylinder(Vec3f center, float radius, float height, int rows, int cols);
		/* Add sphere in space */
		void addSphere(Vec3f center, float radius, int segments, int slices);

		/* Add a simulator to the patches */
		void addNormalSinSimulation();
		/* TODO: Not working yet. */
		void addRandomSphereSimulation();
		/* Add a simulator that rotates the local surfaces */
		void addNormalRotationSimulation();
		void addXYScalingSimulation();
		/* Removes a simulator */
		void removeSimulator(SimulatorTypes simulatorType);
		/* Used for simulation. */
		void update(double dt);

		/* Finalize lattice creation */
		void induceLattice();

		/* Set the name of the lattice */
		void setName(std::string name) { m_name = name; }
		/* Sets wheter the lattice is displayed */
		void setDraw(bool draw) { m_draw = draw; }
		/* Sets wheter the lattice should be animated */
		void setAnimate(bool animate) { m_animate = animate; }
		/* Sets wheter the grid is displayed */
		void setDrawLatticeGrid(bool drawLatticeGrid) { m_drawLatticeGrid = drawLatticeGrid; }
		/* Sets wheter local surfaces are displayed */
		void setDrawLocalSurfaces(bool drawLocalSurfaces) { m_drawLocalSurfaces = drawLocalSurfaces; }
		/* Sets wheter the Blending surface is displayed */
		void setDrawSurface(bool drawSurface) { m_drawSurface = drawSurface; }
		/* Sets wheter normals are displayed */
		bool setDrawNormals(bool drawNormals) { m_drawNormals = drawNormals; }
		/* Sets wireframe mode */
		void setWireframe(bool wireframe) { m_wireframe = wireframe; }
		/* Sets wheter the surface is drawn pixel-accurate */
		void setDrawPixelAccurate(bool drawPixelAccurate) { m_drawPixelAccurate = drawPixelAccurate; }
		/* Set tessellation factors */
		void setTessellationFactors(int inner, int outer) { 
			m_uniforms.tessInner = inner; m_uniforms.tessOuter = outer; }
		/* Set the b function to be used */
		void setBFunction(BFunctionType bFunction) { 
			m_uniforms.bFunctionIndex = static_cast<int>(bFunction); }
		/* Set the matrix */
		void setMatrix(glm::mat4 matrix) { m_matrix = matrix; }
		/* Use different colors per patch */
		void setUseRandomPatchColors(bool useRandomPatchColors) { m_useRandomPatchColors = useRandomPatchColors; }
		/* Set color to be used for every patch, if perPatchColors = false */
		void setPatchColor(glm::vec3 color) { m_color = color; }

		// Returns the name of the Lattice
		std::string name() { return m_name; }

	protected:
		// Called by the update function, can be used in implementations
		virtual void localUpdate(double dt) = 0;

		// Resets all the matrices to their initial matrices. Used to remove any simulation.
		void resetMatrices();

		std::string m_name;
		// The transformation matrix for the Lattice
		glm::mat4 m_matrix;
		glm::mat4 m_view;

		glm::vec3 m_color;
		bool m_useRandomPatchColors = false;

		// Setting m_draw to false should prevent any rendering.
		bool m_draw = true;
		bool m_animate = false;
		bool m_drawLatticeGrid = true;
		bool m_drawLocalSurfaces = false;
		bool m_drawSurface = true;
		bool m_drawNormals = false;
		bool m_wireframe = false;
		bool m_drawPixelAccurate = false;

		std::vector<Locus> m_loci;
		std::vector<Patch> m_patches;

		// Vector of all the control points used by the local surfaces.
		std::vector<glm::vec4> m_controlPoints;
		// TODO: Fix so multiple patches can reuse the same boundaries. To save space.
		std::vector<BoundaryInfo> m_boundaries;
		boundary_map m_boundaryMap;
		size_t m_numUniqueBoundaries = 0;
		// The initial matrices the local surfaces were created with, used for resetting transformations
		std::vector<glm::mat4> m_initialMatrices;
		// The transofrmation matrices for the local surfaces.
		std::vector<glm::mat4> m_matrices;

		// Structure of uniforms used by the shaders.
		struct {
			int tessInner = 10;
			int tessOuter = 10;
			int bFunctionIndex = 0;
			alignas(16) glm::mat4 projection = glm::mat4(1.0f);
			alignas(16) glm::mat4 modelview = glm::mat4(1.0f);
			alignas(16) glm::mat4 normal = glm::mat4(1.0f);
		} m_uniforms;

		int m_simulatorIndex = 0;
		std::unordered_map<SimulatorTypes, std::unordered_map<uint32_t, std::shared_ptr<Simulator>>> m_simulators;
		float m_minAmp = 1.0;
		float m_maxAmp = 10.0;
		float m_minSpeed = 1.0;
		float m_maxSpeed = 5.0;
		float m_minAngle = 10.0;
		float m_maxAngle = 60.0;
		float m_minScale = 0.1f;
		float m_maxScale = 0.5f;

	private:
		// Holds the number of unique points added.
		size_t m_numUniquePointsAdded = 0;
		// Map using the points' pos as a key with a value of the index of the vertex handle of the given point.
		vec3_map m_uniquePointsIndexMap;


		// Check if two points are equal with a tolerance of 1e-5.
		inline bool equal(Vec3f& p1, Vec3f& p2) {
			return std::abs(p1[0] - p2[0]) < 1e-5 &&
				   std::abs(p1[1] - p2[1]) < 1e-5 &&
				   std::abs(p1[2] - p2[2]) < 1e-5;
		}

		/* 
			Used to check if p1 is smaller that p2 along a given direction.
			!!! This might not be correct !!!
		*/
		inline bool smaller(Vec3f& p1, Vec3f& p2, Vec3f& dir) {
			return ((p1[0] - p2[0] * dir[0]) + (p1[1] - p2[1] * dir[1]) + (p1[2] - p2[2] * dir[2])) < 0;
		}

		/*
			Check if face 1 is smaller than face 2 along a given direction
			!!! This is probably not correct !!!
		*/
		inline bool smaller(OpenMesh::FaceHandle& f1, OpenMesh::FaceHandle& f2, Vec3f& dir) {
			Vec3f p1 = point(to_vertex_handle(halfedge_handle(f1)));
			Vec3f p2 = point(to_vertex_handle(halfedge_handle(f2)));
			return smaller(p1, p2, dir);
		}

		// Goes through all the edges and sets the color depending on if it is on the boundary or not.
		void setupEdgeColor();
		// Goes through all vertices and sets the valence, and colors it based on the valence.
		void setupLociValenceAndPointColor();
		// Finds and handles T-loci. Corrects faces and creates local surfaces.
		void handleTLoci();
		// Iterate through vertices and create local surfaces, also sets up patches.
		void setupLocalSurfacesAndPatches();

		void addLatticeProperties();

		void addLocusOnVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
			Vec3f L2RDir, Vec3f T2BDir);
		void addLocusOnCornerVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
			Vec3f L2RDir, Vec3f T2BDir);
		void addLocusOnBoundaryVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
			Vec3f L2RDir, Vec3f T2BDir);
		void addLocusOnInnerVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace,
			Vec3f L2RDir, Vec3f T2BDir);

		void addLocus(OpenMesh::VertexHandle vertex, uint32_t controlPointIndex, uint32_t controlPointCount,
			std::unordered_map<OpenMesh::FaceHandle, uint32_t>& boundaryIndices, Vec3f offset, bool addMatrix = true);

		size_t addBoundaryInfo(BoundaryInfo boundary);

		uint32_t createLocalSurfaceControlPoints(
			Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight);
		uint32_t createLocalSurfaceControlPoints(
			Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
			Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
			Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight);

		// New stuff TODO: ?
		void addLocalSurfaceOnLoci(OpenMesh::VertexHandle vh, Vec3f L2RDir, Vec3f T2BDir);
		void getCornerPointsOfFaceL2RT2B(OpenMesh::FaceHandle& fh, std::vector<Vec3f>& points);

		void addLocalSurfaceOnCornerLocus(OpenMesh::VertexHandle vh, OpenMesh::FaceHandle fh);
		void addLocalSurfaceOnBoundaryLocus(OpenMesh::VertexHandle& vh, std::vector<OpenMesh::FaceHandle>& faces, Vec3f L2RDir, Vec3f T2BDir);
		void addLocalSurfaceOnInnerLocus(OpenMesh::VertexHandle& vh, std::vector<OpenMesh::FaceHandle>& faces, Vec3f L2RDir, Vec3f T2BDir) {}

		// Something random..
		Random m_rng;
	};

}