#ifndef MANIFOLD_MESH_H__
#define MANIFOLD_MESH_H__

#include <vector>
#include <string>

#include <GEL/CGLA/Vec3f.h>

using namespace CGLA;

/** Half-edge data structure for storing a manifold mesh */
class ManifoldMesh {
public:
	typedef int VertKey;
	typedef int FaceKey;
	typedef int EdgeKey;

	static const VertKey INVALID_VERT = -1;
	static const FaceKey INVALID_FACE = -1;
	static const EdgeKey INVALID_EDGE = -1;

	struct Vertex {
		VertKey self;
		EdgeKey edge; // Edge pointing out from this vertex
		Vec3f pos;
		Vec3f normal;

		Vertex() :
			self(INVALID_VERT),
			edge(INVALID_EDGE),
			pos(Vec3f(0)),
			normal(Vec3f(0)) {}
	};

	struct Face {
		FaceKey self;
		EdgeKey edge;
		Vec3f normal;

		Face() :
			self(INVALID_FACE),
			edge(INVALID_EDGE),
			normal(Vec3f(0)) {}
	};

	struct HalfEdge {
		EdgeKey self;
		EdgeKey prev;
		EdgeKey next;
		EdgeKey twin;

		VertKey vert; // Vertex this edge points out of
		FaceKey face;

		HalfEdge() :
			self(INVALID_EDGE),
			prev(INVALID_EDGE),
			next(INVALID_EDGE),
			twin(INVALID_EDGE),
			vert(INVALID_VERT),
			face(INVALID_FACE) {}
	};

	std::vector<Vertex> vertices;
	std::vector<Face> faces;
	std::vector<HalfEdge> edges;

	explicit ManifoldMesh() = default;
	ManifoldMesh(const ManifoldMesh& other);
	ManifoldMesh(ManifoldMesh&& other);

	ManifoldMesh& operator=(const ManifoldMesh& other);
	ManifoldMesh& operator=(ManifoldMesh&& other);

	void clear();
	void build(std::vector<Vec3f> vertexPositions,
		const std::vector<std::vector<int>>& faceIdxLists,
		size_t expectedEdgeCount = 0);

	void laplaceSmooth(float lambda, int niter = 1);
	void taubinSmooth(float lambda, float mu, int niter = 1);

	void splitEdge(EdgeKey ek, float a = 0.5f);
	void flipEdge(EdgeKey ek);
	void removeEdge(EdgeKey ek);

	void removeInvalid();

	float computeBoundingSphere() const noexcept;
	float computeBoundingSphere(Vec3f& center) const noexcept;

	std::vector<float> gaussCurvatures() const;

	Vec3f faceNormal(FaceKey fk) const;
	Vec3f faceNormal(const Face& f) const;
	void computeFaceNormals(bool normalize = true);

	Vec3f vertexNormal(VertKey vk, bool faceNormalsComputed = false) const;
	Vec3f vertexNormal(const Vertex& v, bool faceNormalsComputed = false) const;
	void computeVertexNormals(bool faceNormalsComputed = false, bool normalize = true);

	void saveToObj(const std::string& fname) const;
	void loadFromObj(const std::string& fname);

	EdgeKey next(EdgeKey ek) const;
	EdgeKey prev(EdgeKey ek) const;
	EdgeKey twin(EdgeKey ek) const;

	Vec3f& vpos(const Vertex& v);
	const Vec3f& vpos(const Vertex& v) const;
	Vec3f& vpos(VertKey vk);
	const Vec3f& vpos(VertKey vk) const;

	void assertConsistent() const;
};

bool operator==(const ManifoldMesh& lhs, const ManifoldMesh& rhs);
bool operator!=(const ManifoldMesh& lhs, const ManifoldMesh& rhs);

bool operator==(const ManifoldMesh::Vertex& lhs, const ManifoldMesh::Vertex& rhs);
bool operator!=(const ManifoldMesh::Vertex& lhs, const ManifoldMesh::Vertex& rhs);

bool operator==(const ManifoldMesh::Face& lhs, const ManifoldMesh::Face& rhs);
bool operator!=(const ManifoldMesh::Face& lhs, const ManifoldMesh::Face& rhs);

bool operator==(const ManifoldMesh::HalfEdge& lhs, const ManifoldMesh::HalfEdge& rhs);
bool operator!=(const ManifoldMesh::HalfEdge& lhs, const ManifoldMesh::HalfEdge& rhs);

#endif // MANIFOLD_MESH_H__
