#pragma once

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <unordered_map>

namespace OML {

	using Vec2f = OpenMesh::Vec2f;
	using Vec3f = OpenMesh::Vec3f;
	using Col3 = OpenMesh::Vec3uc;

	struct BoundaryInfo
	{
		BoundaryInfo()
			: us(0.0f), ue(1.0f), vs(0.0f), ve(1.0f) {}
		BoundaryInfo(float us, float ue, float vs, float ve)
			: us(us), ue(ue), vs(vs), ve(ve) {}

		float us, ue, vs, ve;
	};

	enum class LocalSurfaceType
	{
		Bezier3x3 = 0
	};

	struct Locus
	{
		Locus() : controlPointIndex(0), controlPointCount(0), matrixIndex(0) {}

		uint32_t controlPointIndex;
		uint32_t controlPointCount;
		uint32_t matrixIndex;
		OpenMesh::VertexHandle vh;
		std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo> faceMappings;
	};

	struct Patch
	{
		std::array<size_t, 4> lociIndices;
		OpenMesh::FaceHandle fh;
	};

	enum class BFunctionType {
		B1Poly = 0, B2Poly, LERBS
	};

	// Custom traits
	struct LatticeTraits : public OpenMesh::DefaultTraits
	{
		typedef Vec3f Point;

		VertexAttributes(OpenMesh::Attributes::Color | OpenMesh::Attributes::Status); // Use color for vailence

		EdgeAttributes(OpenMesh::Attributes::Color | OpenMesh::Attributes::Status); // Use color for on/inside boundary

		HalfedgeAttributes(OpenMesh::Attributes::Status); // Use for deleting half edges

		FaceAttributes(OpenMesh::Attributes::Status); // Use for deleting faces
	};

	namespace LatticeProperties {
		static OpenMesh::VPropHandleT<size_t> VertexValence;
		static OpenMesh::VPropHandleT<uint32_t> LocusIndex;
	};

	const Col3 BOUNDARY_EDGE_COLOR = Col3(0, 0, 0);
	const Col3 INNER_EDGE_COLOR = Col3(200, 200, 200);
	const std::vector<Col3> LOCUS_VALENCE_COLOR = {
		Col3(0,255,255), Col3(255,255,255), Col3(255,0,0), Col3(0,255,0), Col3(0,0,255), Col3(0,255,255)
	};	

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

	protected:
		virtual void setupLocalSurfaceVertex(Locus& locus) = 0;
		virtual void setupPatchVertices(Patch& patch) = 0;

		// Lattice stuff
		std::string m_name;

		bool m_draw = true;
		bool m_animate = false;
		bool m_drawLatticeGrid = true;
		bool m_drawLocalSurfaces = false;
		bool m_drawSurface = true;
		bool m_drawNormals = false;
		bool m_wireframe = false;
		bool m_drawPixelAccurate = false;

		uint32_t m_numControlPoints;
		uint32_t m_numLoci;
		uint32_t m_numPatches;

		std::vector<Locus> m_loci;
		std::vector<Patch> m_patches;

		std::vector<glm::vec4> m_controlPoints;
		std::vector<BoundaryInfo> m_boundaries;
		std::vector<glm::mat4> m_matrices;

		struct {
			int tessInner = 10;
			int tessOuter = 10;
			int bFunctionIndex = 0;
			alignas(16) glm::mat4 projection = glm::mat4(1.0f);
			alignas(16) glm::mat4 modelview = glm::mat4(1.0f);
		} m_uniforms;

	private:
		// Helper functions for setting up the lattice stuff
		void setupEdgeColor();
		void setupLociValenceAndPointColor();
		void handleTLoci();
		void setupLocalSurfacesAndPatches();

		void addLatticeProperties();

		void addLocusOnVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);
		void addLocusOnCornerVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);
		void addLocusOnBoundaryVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);
		void addLocusOnInnerVertex(
			OpenMesh::VertexHandle vertex, OpenMesh::VertexHandle next_vertex,
			OpenMesh::VertexHandle prev_vertex, int vertexIndexOnFace);

		void addLocus(OpenMesh::VertexHandle vertex, std::vector<Vec3f>& controlPoints,
			std::unordered_map<OpenMesh::FaceHandle, BoundaryInfo>& faceMappings, Vec3f offset);

		std::vector<Vec3f> createLocalSurfaceControlPoints(
			Vec3f topLeft, Vec3f topRight, Vec3f bottomLeft, Vec3f bottomRight);
		std::vector<Vec3f> createLocalSurfaceControlPoints(
			Vec3f topLeft, Vec3f topMiddle, Vec3f topRight,
			Vec3f middleLeft, Vec3f middle, Vec3f middleRight,
			Vec3f bottomLeft, Vec3f bottomMiddle, Vec3f bottomRight);
	};

}