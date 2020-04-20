#include "subdivided_icosahedron.h"
#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <array>

SubdividedIcosahedron::SubdividedIcosahedron(float r) :
    ManifoldMesh(),
    center_(0.0f),
    r_(r),
    divisionLevel_(1)
{
    buildIcosahedron();
}

SubdividedIcosahedron::SubdividedIcosahedron(const Vec3f& center, float r, int divLvl) :
    ManifoldMesh(),
    center_(center),
    r_(r),
    divisionLevel_(divLvl)
{
    buildIcosahedron();
    if (divLvl > 1) {
        subdivide(divLvl);
    }
}

void SubdividedIcosahedron::buildIcosahedron()
{
    buildIcosahedron(center_, r_);
}

void SubdividedIcosahedron::buildIcosahedron(const Vec3f& center, float r)
{
    // Make sure everything is cleared if we already had a mesh
    clear();

    // Set fields
    this->center_ = center;
    this->r_ = r;
    divisionLevel_ = 1;

    // Hard coded coords. for unit icosahedron
    constexpr float c1 = 0.0;
    constexpr float c2 = static_cast<float>(0.52573111211913360602566908484786);
    constexpr float c3 = static_cast<float>(0.85065080835203993218154049706296);

    std::vector<Vec3f> vertexPositions = {
        Vec3f(-c3, c1, -c2),
        Vec3f(-c3, c1, c2),
        Vec3f(-c2, -c3, c1),
        Vec3f(-c2, c3, c1),
        Vec3f(c1, -c2, -c3),
        Vec3f(c1, -c2, c3),
        Vec3f(c1, c2, -c3),
        Vec3f(c1, c2, c3),
        Vec3f(c2, -c3, c1),
        Vec3f(c2, c3, c1),
        Vec3f(c3, c1, -c2),
        Vec3f(c3, c1, c2)
    };

    static const std::vector<std::vector<int>> faceIdxLists = {
        { 0, 1, 3 },
        { 1, 5, 7 },
        { 1, 7, 3 },
        { 4, 10, 8 },
        { 0, 3, 6 },
        { 0, 6, 4 },
        { 4, 6, 10 },
        { 0, 4, 2 },
        { 1, 2, 5 },
        { 0, 2, 1 },
        { 2, 8, 5 },
        { 2, 4, 8 },
        { 5, 8, 11 },
        { 5, 11, 7 },
        { 8, 10, 11 },
        { 3, 9, 6 },
        { 3, 7, 9 },
        { 6, 9, 10 },
        { 9, 11, 10 },
        { 7, 11, 9 }
    };

    // Scale and translate vertices
    for (Vec3f& v : vertexPositions) {
        v *= r;
        v += center;
    }

    // Build mesh
    ManifoldMesh::build(vertexPositions, faceIdxLists, 2*30);
}

void SubdividedIcosahedron::rescale(float newR) noexcept
{
    const float scale = newR / r_;
    for (auto& v : vertices) {
        v.pos = (v.pos - center_)*scale + center_;
    }
    r_ = newR;
}

void SubdividedIcosahedron::move(const Vec3f& newCenter) noexcept
{
    Vec3f diff = newCenter - center_;
    for (auto& v : vertices) {
        v.pos += diff;
    }
    center_ = newCenter;
}

void SubdividedIcosahedron::moveAndRescale(const Vec3f& newCenter, float newR) noexcept
{
    const float scale = newR / r_;
    for (auto& v : vertices) {
        v.pos = (v.pos - center_)*scale + newCenter;
    }
    r_ = newR;
    center_ = newCenter;
}

void SubdividedIcosahedron::subdivide(int numDivide)
{
    for (int i = 0; i < numDivide; i++) {
        singleSubdivide_();
    }
    divisionLevel_ += numDivide;
}

float SubdividedIcosahedron::r() const noexcept
{
    return r_;
}

const Vec3f& SubdividedIcosahedron::center() const noexcept
{
    return center_;
}

void SubdividedIcosahedron::singleSubdivide_()
{
    std::vector<Face> newFaces;

    newFaces.reserve(4 * faces.size());
	vertices.reserve(vertices.size() + edges.size() / 2);
    edges.reserve(2 * edges.size() + 6 * faces.size());

    std::unordered_map<EdgeKey, VertKey> addedVertices;
    addedVertices.reserve(edges.size());
    for (Face& f : faces) {
        // Add new faces
        std::array<FaceKey, 4> extraFaces;
        for (int i = 0; i < 4; i++) {
            extraFaces[i] = newFaces.size();
            newFaces.push_back(Face());
            newFaces.back().self = extraFaces[i];
        }

        std::array<EdgeKey, 3> oldEdges;
        std::array<EdgeKey, 3> newEdges;

        // First, split all old edges
        EdgeKey ek = f.edge;
        for (int i = 0; i < 3; i++) {
            HalfEdge& e = edges[ek];
            Vertex& v = vertices[e.vert];

            // Make new edge
            EdgeKey ekn = edges.size();
            edges.push_back(HalfEdge());
            HalfEdge& en = edges.back();

            oldEdges[i] = ek;
            newEdges[i] = ekn;

            // Get new vertex
            VertKey vkn;
            const auto findRes = addedVertices.find(ek);
            if (findRes == addedVertices.end()) {
                // New vertex not yet created so make it now
                vkn = vertices.size();
                vertices.push_back(Vertex());
                Vertex& vn = vertices.back();
                vn.edge = ekn;
                vn.self = vkn;
                Vertex vnOld = vertices[edges[e.next].vert];
                vn.pos = normalize(0.5*(v.pos + vnOld.pos) - center_)*r_ + center_;

                // We have not processed this edge so the twin pointer for e still points to an old edge
                en.twin = e.twin;
                edges[en.twin].twin = ekn;

                addedVertices[ek] = vkn;
                addedVertices[twin(ek)] = vkn;
            } else {
                // New vertex already made so retrieve its key
                vkn = findRes->second;

                // We have already processed this edge so the twin pointer for e points to a new edge
                en.twin = prev(twin(prev(twin(prev(twin(ek))))));
                edges[en.twin].twin = ekn;
            }

            // Set fields for new edge
            en.self = ekn;
            en.next = e.next;
            en.vert = vkn;
            en.face = extraFaces[(i + 1) % 3];

            edges[en.next].prev = ekn;
            e.face = extraFaces[i % 3];

            // Register edge with new face
            newFaces[extraFaces[i % 3]].edge = ek;

            // Prepare for next iteration
            ek = en.next;
        }

        // Add half-edges for outermost new faces
        for (int i = 0; i < 3; i++) {
            // Make new half-edge
            EdgeKey ekn = edges.size();
            edges.push_back(HalfEdge());
            HalfEdge& en = edges.back();

            // Set fields - twin is set later
            en.self = ekn;
            en.prev = oldEdges[i];
            en.next = newEdges[(i + 2) % 3];
            en.vert = edges[newEdges[i]].vert;
            en.face = extraFaces[i];
            edges[oldEdges[i]].next = ekn;
            edges[newEdges[(i + 2) % 3]].prev = ekn;
        }

        // Add half-edges for remaning innermost face
        std::array<EdgeKey, 3> innerEdges;
        for (int i = 0; i < 3; i++) {
            innerEdges[i] = edges.size();
            edges.push_back(HalfEdge());
        }
        for (int i = 0; i < 3; i++) {
            HalfEdge& e = edges[innerEdges[i]];
            e.self = innerEdges[i];
            e.prev = innerEdges[(i + 2) % 3];
            e.next = innerEdges[(i + 1) % 3];
            e.twin = next(oldEdges[i]);
            e.vert = edges[prev(oldEdges[i])].vert;
            e.face = extraFaces[3];

            edges[next(oldEdges[i])].twin = innerEdges[i];
        }
        newFaces[extraFaces[3]].edge = innerEdges[0];
    }
    faces = std::move(newFaces);
}