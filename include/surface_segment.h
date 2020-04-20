#ifndef SURFACE_SEGMENT_H__
#define SURFACE_SEGMENT_H__

#include <vector>
#include <unordered_set>
#include <utility>

#include "graph.h"
#include <GEL/CGLA/Vec3f.h>

#include "manifold_mesh.h"
#include "subdivided_icosahedron.h"
#include "volume.h"

using namespace CGLA;

using FloatGraph = Graph<float, float, float>;

enum CostType : int {
	ON_SURFACE = 0,
	IN_SURFACE = 1
};

ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType);
ManifoldMesh surfaceCut(const Volume<float>& cost, const ManifoldMesh&& init,
    int numSamples, float sampleStep, int maxDiff, CostType costType);

std::pair<std::vector<ManifoldMesh>, size_t> surfaceCutPlaneSepQPBO(const Volume<float>& cost,
	std::vector<ManifoldMesh> meshes, int numSamples, float sampleStep, int maxDiff, CostType costType,
	const std::vector<Vec3f>& centers, const std::vector<std::vector<size_t>>& connections);

#endif // SURFACE_SEGMENT_H__