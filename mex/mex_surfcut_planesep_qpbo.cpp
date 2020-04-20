#include <vector>
#include <thread>
#include <memory>
#include <stdexcept>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/Vec3f.h>

#include "util.h"
#include "volume.h"
#include "subdivided_icosahedron.h"
#include "surface_segment.h"
#include "matlab_util.h"

using namespace CGLA;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ensureOrError(9 <= nrhs && nrhs <= 11, "Must supply between 7 and 11 inputs");
	Volume<float> cost = getVolumeChecked<float>(prhs[0], "Cost volume");
	Volume<float> centers = getCastVolumeChecked<float>(prhs[1], "Centers");
	std::shared_ptr<float> initRs = getCastSharedPtrChecked<float>(prhs[2], "initR");
	int initSubDiv = getCastScalarChecked<int>(prhs[3], "initSubDiv");
	int numSamples = getCastScalarChecked<int>(prhs[4], "numSamples");
	float sampleStep = getCastScalarChecked<float>(prhs[5], "sampleStep");
	int maxDiff = getCastScalarChecked<int>(prhs[6], "maxDiff");
	const mxArray *connections = prhs[7];
	CostType costType = static_cast<CostType>(getCastScalarChecked<int>(prhs[8], "costType"));
	float smoothFactor = nrhs > 9 ? getCastScalarChecked<float>(prhs[9], "smoothFactor") : 0.0f;
	int smoothIter = nrhs > 10 ? getCastScalarChecked<int>(prhs[10], "smoothIter") : 0;

	// TODO: Extra input validation

	ensureOrError(centers.ny == 3 && isMatrix(prhs[1]), "Centers must be N x 3");
	ensureOrError(mxIsScalar(prhs[2]) || isVector(prhs[2]), "initR must be a scalar or vector");
	ensureOrError(mxIsSparse(connections), "Connections must be a sparse matrix");
	ensureOrError(!mxIsComplex(connections), "Connections must be real");
	ensureOrError(mxGetN(connections) == mxGetM(connections), "Connections must be a square matrix");

	int nr = mxGetNumberOfElements(prhs[2]);
	int nmesh = centers.nx;
	ensureOrError(nr == 1 || nr == nmesh, "Must provide a single initR or one for every center");
	ensureOrError(mxGetN(connections) == nmesh, "Connections must have a row for each mesh");

	std::vector<ManifoldMesh> meshes;
	std::vector<float *> vertData;
	std::vector<int *> faceData;
	std::vector<std::vector<size_t>> connVec;
	std::vector<Vec3f> centerVec;

	// Allocate all needed memory and prepare for cut
	mxArray *vertCell = mxCreateCellMatrix(1, nmesh);
	mxArray *faceCell = mxCreateCellMatrix(1, nmesh);
	meshes.reserve(nmesh);
	vertData.reserve(nmesh);
	faceData.reserve(nmesh);
	connVec.reserve(nmesh);
	float *initR = initRs.get();
	for (int i = 0; i < nmesh; ++i) {
		Vec3f center(centers.at(i, 0, 0), centers.at(i, 1, 0), centers.at(i, 2, 0));
		meshes.push_back(SubdividedIcosahedron(center, *initR, initSubDiv));

		const size_t numVerts = meshes[i].vertices.size();
		const size_t numFaces = meshes[i].faces.size();
		mxArray *vertices = mxCreateNumericMatrix(numVerts, 3, mxSINGLE_CLASS, mxREAL);
		mxArray *faces = mxCreateNumericMatrix(numFaces, 3, mxINT32_CLASS, mxREAL);
		vertData.push_back(static_cast<float *>(mxGetData(vertices)));
		faceData.push_back(static_cast<int *>(mxGetData(faces)));

		mxSetCell(vertCell, i, vertices);
		mxSetCell(faceCell, i, faces);

		centerVec.push_back(center);

		// When I wrote this, only MathWorks and I knew how this worked.
		// Now, only MathWorks knows...
		size_t *ir = mxGetIr(connections);
		size_t *jc = mxGetJc(connections);
		size_t nconn = jc[i + 1] - jc[i];
		std::vector<size_t> conni;
		conni.reserve(nconn);
		for (size_t ci = jc[i]; ci < jc[i + 1]; ++ci) {
			conni.push_back(ir[ci]);
		}
		connVec.push_back(conni);

		if (nr > 1) {
			initR++;
		}
	}

	// Do the cut
	size_t totalUnlabeled;
	std::tie(meshes, totalUnlabeled) = surfaceCutPlaneSepQPBO(cost, std::move(meshes),
		numSamples, sampleStep, maxDiff, costType, centerVec, connVec);
	/*try {
	} catch (const std::runtime_error& err) {
		mexErrMsgIdAndTxt("surfseg:internal", err.what());
	}*/
	// Run smoothing
	for (auto& m : meshes) {
		m.taubinSmooth(smoothFactor, 1.04f * smoothFactor, smoothIter);
	}
	// Extract the data
	for (int i = 0; i < nmesh; ++i) {
		const auto& m = meshes[i];
		const size_t numVerts = meshes[i].vertices.size();
		const size_t numFaces = meshes[i].faces.size();

		for (const auto& v : m.vertices) {
			vertData[i][v.self + 0 * numVerts] = v.pos[0];
			vertData[i][v.self + 1 * numVerts] = v.pos[1];
			vertData[i][v.self + 2 * numVerts] = v.pos[2];
		}

		for (const auto& f : m.faces) {
			// Add one since MATLAB uses 1-indexing
			faceData[i][f.self + 0 * numFaces] = m.edges[f.edge].vert + 1;
			faceData[i][f.self + 1 * numFaces] = m.edges[m.next(f.edge)].vert + 1;
			faceData[i][f.self + 2 * numFaces] = m.edges[m.next(m.next(f.edge))].vert + 1;
		}
	}

	plhs[0] = faceCell;
	plhs[1] = vertCell;
	plhs[2] = mxCreateDoubleScalar(totalUnlabeled);
}
