#include "surface_segment.h"
#include "util.h"
#include "QPBO.h"
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <sstream>


ManifoldMesh& updateVertices(const FloatGraph& graph, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t k = 0, size_t offset = 0);

size_t updateVerticesQPBO(const QPBO<float>& qpbo, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t k = 0, size_t offset = 0);

FloatGraph& buildSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k = 0, size_t offset = 0);

QPBO<float>& buildQPBOSurfaceGraph(QPBO<float>& qpbo, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k = 0, size_t offset = 0);


Volume<float>& extractCostSamples(const Volume<float>& cost, const ManifoldMesh& mesh,
	Volume<float>& samples, int numSamples, float sampleStep, size_t k = 0);

ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType)
{
    ManifoldMesh mesh(init);
    return surfaceCut(cost, std::move(mesh), numSamples, sampleStep, maxDiff, costType);
}

ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh&& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType)
{
    ManifoldMesh mesh(std::move(init));
    mesh.computeVertexNormals(false, true); // Ensure these are correct
    size_t numNodes = init.vertices.size();

    // Make cost sample volume
    Volume<float> costSamples(numNodes, numSamples, 1);
    costSamples.alloc();
    extractCostSamples(cost, mesh, costSamples, numSamples, sampleStep);

    // Build min-cut graph and find optimal cut
    size_t totalSamples = costSamples.numElem();
    FloatGraph graph(totalSamples,
        mesh.edges.size() * (numSamples - maxDiff - 1) + totalSamples - costSamples.nx,
		graphErrFunc
	);
    graph.add_node(totalSamples);
    buildSurfaceGraph(graph, costSamples, mesh, maxDiff, costType);

    graph.maxflow();

    // Update mesh vertex positions
	updateVertices(graph, costSamples, mesh, sampleStep);

    return mesh;
}

std::pair<std::vector<ManifoldMesh>, size_t> surfaceCutPlaneSepQPBO(const Volume<float>& cost,
	std::vector<ManifoldMesh> meshes, int numSamples, float sampleStep, int maxDiff, CostType costType,
	const std::vector<Vec3f>& centers, const std::vector<std::vector<size_t>>& connections)
{
	assert(meshes.size() == centers.size());
	assert(meshes.size() == connections.size());
	for (auto& m : meshes) {
		m.computeVertexNormals(false, true); // Ensure these are correct
	}

	// Make cost sample volumes
	std::vector<Volume<float>> costSamples;
	costSamples.reserve(meshes.size());
	size_t totalSamples = 0;
	size_t totalEdges = 0;
	for (const auto& m : meshes) {
		size_t numNodes = m.vertices.size();
		Volume<float> cs(numNodes, numSamples, 1);
		cs.alloc();
		extractCostSamples(cost, m, cs, numSamples, sampleStep);
		size_t layerLen = numNodes * numSamples;

		totalSamples += cs.numElem();
		totalEdges += m.edges.size() * (numSamples - maxDiff - 1) + cs.numElem() - cs.nx;
		costSamples.push_back(std::move(cs));
	}

	// Compute number of needed nodes for the surface to plane stuff
	size_t extraNodes = 0;
	for (int i = 0; i < meshes.size(); ++i) {
		Vec3f ceni = centers[i];
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			extraNodes += roundf(length(ceni - cenj) / sampleStep);
		}
	}

	// Build surface graphs
	std::vector<size_t> offsets;
	offsets.reserve(meshes.size());
	// Adding an extra totalEdges is a heuristic for how many surface to plane position edges we are gonna need
	// It seems to be a good upper bound, but there are no guarantees
	QPBO<float> qpbo(totalSamples + extraNodes, 2 * totalEdges, graphErrFunc);;

	//graph.add_node(totalSamples + extraNodes);
	qpbo.AddNode(totalSamples + extraNodes);
	for (size_t i = 0, offset = 0; i < meshes.size(); ++i) {
		auto& cs = costSamples[i];
		const auto& m = meshes[i];

		buildQPBOSurfaceGraph(qpbo, cs, m, maxDiff, costType, 0, offset);

		offsets.push_back(offset);
		offset += cs.numElem();
		// From this point we don't need the actual samples anymore, but we want to keep the Volume
		// for computing indices
		cs.data = nullptr;
	}

	// Add plane edges
	int totalPlaneArcs = 0;
	for (int i = 0, planeNodeOffset = totalSamples; i < meshes.size(); ++i) {
		Vec3f ceni = centers[i];
		for (int j : connections[i]) {
			Vec3f cenj = centers[j];
			int numPlaneNodes = roundf(length(ceni - cenj) / sampleStep);
			Vec3f normal = normalize(cenj - ceni);
			// Add intracolumn plane edges
			for (int n = 0; n < numPlaneNodes; ++n) {
				if (n > 0) {
					size_t nip = planeNodeOffset + n;
					size_t njp = planeNodeOffset + n - 1;
					qpbo.AddPairwiseTerm(nip, njp, 0, infOrMax<float>(), 0, 0);
					totalPlaneArcs++;
				}
			}
			// Add edges from surface to plane
			for (const auto& v : meshes[i].vertices) {
				Vec3f p = v.pos;
				/*if (dot(normal, p) <= dot(normal, ceni) && dot(normal, v.normal) <= 0) {
					// Vertex is behind the plane and normal points away, so we can skip it
					continue;
				}*/
				for (int n = 0, si = 0; n < numPlaneNodes; ++n) {
					// We don't reset si, since if a point was behind the plane for some n, then it will also
					// be behind it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) >= dist) {
							size_t nip = costSamples[i].idx(v.self, si, 0) + offsets[i];
							size_t njp = planeNodeOffset + n;
							qpbo.AddPairwiseTerm(nip, njp, 0, infOrMax<float>(), 0, 0);
							totalPlaneArcs++;
							break;
						}
						p += sampleStep * v.normal;
					}
				}
			}
			// Add edges from plane to surface
			for (const auto& v : meshes[j].vertices) {
				Vec3f p = v.pos;
				/*if (dot(normal, p) >= dot(normal, ceni) && dot(normal, v.normal) >= 0) {
					// Vertex is in front of the plane and normal points away, so we can skip it
					continue;
				}*/
				for (int n = numPlaneNodes - 1, si = 0; n >= 0; --n) {
					// We don't reset si, since if a point was in front of the plane for some n, then it will
					// also be in front of it for n + 1
					float dist = dot(normal, ceni + n * sampleStep * normal);
					for (; si < numSamples; ++si) {
						if (dot(normal, p) <= dist) {
							size_t nip = costSamples[j].idx(v.self, si, 0) + offsets[j];
							size_t njp = planeNodeOffset + n;
							qpbo.AddPairwiseTerm(nip, njp, infOrMax<float>(), 0, 0, 0);
							totalPlaneArcs++;
							break;
						}
						p += sampleStep * v.normal;
					}
				}
			}
			planeNodeOffset += numPlaneNodes;
		}
	}

	qpbo.Solve();
	qpbo.ComputeWeakPersistencies();

	// Update mesh vertex positions
	size_t totalUnlabeled = 0;
	for (size_t i = 0; i < meshes.size(); ++i) {
		const auto& cs = costSamples[i];
		auto& m = meshes[i];
		totalUnlabeled += updateVerticesQPBO(qpbo, cs, m, sampleStep, 0, offsets[i]);
	}
	return std::make_pair(meshes, totalUnlabeled);
}

ManifoldMesh& updateVertices(const FloatGraph& graph, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t k, size_t offset)
{
	for (auto& v : mesh.vertices) {
		// Find upper position for this vertex
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			if (graph.what_segment(ni) == SOURCE) {
				//v.pos += i * sampleStep * v.normal * Vec3f(1, 1, 0.09); // DEBUG
				//v.pos += i * sampleStep * v.normal * Vec3f(1, 1, 0.625); // DEBUG
				v.pos += i * sampleStep * v.normal;
				break;
			}
		}
	}

	return mesh;
}

size_t updateVerticesQPBO(const QPBO<float>& qpbo, const Volume<float>& costSamples,
	ManifoldMesh& mesh, float sampleStep, size_t k, size_t offset)
{
	size_t numUnlabeled = 0;
	for (auto& v : mesh.vertices) {
		// Find upper position for this vertexfor (int i = costSamples.ny - 1; i >= 0; --i) {
		for (int i = costSamples.ny - 1; i >= 0; --i) {
			size_t ni = costSamples.idx(v.self, i, k) + offset;
			int label = qpbo.GetLabel(ni);
			if (label < 0) {
				numUnlabeled++;
			} else if (label == 0) {
				//v.pos += i * sampleStep * v.normal * Vec3f(1, 1, 0.09); // DEBUG
				//v.pos += i * sampleStep * v.normal * Vec3f(1, 1, 0.625); // DEBUG
				v.pos += i * sampleStep * v.normal;
				break;
			}
		}
	}
	/*if (numUnlabeled > 0) {
		// Some nodes were unlabeled so abort with an error
		std::ostringstream oss;
		oss << numUnlabeled << "/" << qpbo.GetNodeNum();
		oss << " (" << 100.0 * static_cast<double>(numUnlabeled) / qpbo.GetNodeNum() << "%)";
		oss << " nodes were unlabeled";
		throw std::runtime_error(oss.str());
	}*/

	return numUnlabeled;
}

FloatGraph& buildSurfaceGraph(FloatGraph& graph, const Volume<float>& costSamples,
    const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k, size_t offset)
{
    using VertKey = ManifoldMesh::VertKey;
    using EdgeKey = ManifoldMesh::EdgeKey;

    size_t numSamples = costSamples.ny;

    // Add edges
    for (const auto& v : mesh.vertices) {
        // Add intracolumn (downward) edges
        for (int i = numSamples - 1; i > 0; --i) {
            size_t ni = costSamples.idx(v.self, i, k) + offset;
            size_t nj = costSamples.idx(v.self, i - 1, k) + offset;
            graph.add_edge(ni, nj, infOrMax<float>(), 0);
        }

        // Add intercolumn (neighbor) edges
        EdgeKey ek0 = mesh.twin(v.edge);
        EdgeKey ek = ek0;
        do {
            assert(mesh.edges[ek].vert != v.self && mesh.edges[mesh.twin(ek)].vert == v.self);
            VertKey nk = mesh.edges[ek].vert;
            for (int i = numSamples - 1; i > maxDiff; --i) {
                size_t ni = costSamples.idx(v.self, i, k) + offset;
                size_t nj = costSamples.idx(nk, i - maxDiff, k) + offset;
                graph.add_edge(ni, nj, infOrMax<float>(), 0);
            }
            ek = mesh.twin(mesh.next(ek)); // Next neighbor
        } while (ek != ek0);

        // Add source and sink edges
        float prevC = 0.0f;
        for (int i = 0; i < numSamples; ++i) {
            size_t ni = costSamples.idx(v.self, i, k) + offset;
            float c = costSamples.at(v.self, i, k);
            float w = costType == ON_SURFACE ? c - prevC : c;
			if (i == 0) {
				w = -1.0f;
			}
            if (w < 0) {
                graph.add_tweights(ni, -w, 0);
            } else {
                graph.add_tweights(ni, 0, w);
            }
            prevC = c;
        }
    }

    return graph;
}

QPBO<float>& buildQPBOSurfaceGraph(QPBO<float>& qpbo, const Volume<float>& costSamples,
	const ManifoldMesh& mesh, int maxDiff, CostType costType, size_t k, size_t offset)
{
	// TODO: Merge with buildSurfaceGraph
	using VertKey = ManifoldMesh::VertKey;
	using EdgeKey = ManifoldMesh::EdgeKey;

	size_t numSamples = costSamples.ny;
	size_t totalSamples = costSamples.numElem();

	// Add edges
	for (const auto& v : mesh.vertices) {
		// Add intracolumn (downward/upward) edges
		for (int i = numSamples - 1; i > 0; --i) {
			size_t nip = costSamples.idx(v.self, i, k) + offset;
			size_t njp = costSamples.idx(v.self, i - 1, k) + offset;
			qpbo.AddPairwiseTerm(nip, njp, 0, infOrMax<float>(), 0, 0);
		}

		// Add intercolumn (neighbor) edges
		EdgeKey ek0 = mesh.twin(v.edge);
		EdgeKey ek = ek0;
		do {
			assert(mesh.edges[ek].vert != v.self && mesh.edges[mesh.twin(ek)].vert == v.self);
			VertKey nk = mesh.edges[ek].vert;
			for (int i = numSamples - 1; i > maxDiff; --i) {
				size_t nip = costSamples.idx(v.self, i, k) + offset;
				size_t njp = costSamples.idx(nk, i - maxDiff, k) + offset;
				qpbo.AddPairwiseTerm(nip, njp, 0, infOrMax<float>(), 0, 0);
			}
			ek = mesh.twin(mesh.next(ek)); // Next neighbor
		} while (ek != ek0);

		// Add source and sink edges
		float prevC = 0.0f;
		for (int i = 0; i < numSamples; ++i) {
			size_t nip = costSamples.idx(v.self, i, k) + offset;
			float c = costSamples.at(v.self, i, k);
			float w = (costType == ON_SURFACE ? c - prevC : c);
			if (i == 0) {
				w = -1.0f;
			}
			if (w < 0) {
				qpbo.AddUnaryTerm(nip, 0, -w);
			} else {
				qpbo.AddUnaryTerm(nip, w, 0);
			}
			prevC = c;
		}
	}

	return qpbo;
}

Volume<float>& extractCostSamples(const Volume<float>& cost, const ManifoldMesh& mesh,
    Volume<float>& samples, int numSamples, float sampleStep, size_t k)
{
    // Sample volume along vertex normals
    for (const auto& v : mesh.vertices) {
        Vec3f p = v.pos;
        for (int i = 0; i < numSamples; ++i) {
			//p[2] = fmaxf(0.0f, fminf(cost.nz - 1, p[2])); // Pad volume to repeat top and bottom slice
			//samples.at(v.self, i, k) = cost.contains(p) ? cost.interp(p) : infOrMax<float>(); // Extrapolate inf
			samples.at(v.self, i, k) = cost.contains(p) ? cost.interp(p) : 0; // Extrapolate 0
            //p += sampleStep * v.normal * Vec3f(1, 1, 0.09); // DEBUG
			//p += sampleStep * v.normal * Vec3f(1, 1, 0.625); // DEBUG
			p += sampleStep * v.normal;
        }
    }

    return samples;
}
