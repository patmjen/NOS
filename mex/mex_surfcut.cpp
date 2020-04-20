#include <vector>
#include <thread>
#include <memory>
#include <cmath>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Mat3x3f.h>

#include "volume.h"
#include "surface_segment.h"
#include "matlab_util.h"

using namespace CGLA;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	ensureArgRange(nrhs, 8, 10);
    ensureSize(prhs[1], { -1, 3 }, "Centers", "N");
    ensureOrError(mxIsScalar(prhs[2]) || isVector(prhs[2]) || isSize(prhs[2], { 3, 3, -1 }),
        "initR must be a scalar or vector or 3 x 3 x N");
	Volume<float> cost = getVolumeChecked<float>(prhs[0], "Cost volume");
	Volume<float> centers = getCastVolumeChecked<float>(prhs[1], "Centers");
	Volume<float> initRs = getCastVolumeChecked<float>(prhs[2], "initRs");
	int initSubDiv = getCastScalarChecked<int>(prhs[3], "initSubDiv");
	int numSamples = getCastScalarChecked<int>(prhs[4], "numSamples");
	float sampleStep = getCastScalarChecked<float>(prhs[5], "sampleStep");
	int maxDiff = getCastScalarChecked<int>(prhs[6], "maxDiff");
	CostType costType = static_cast<CostType>(getCastScalarChecked<int>(prhs[7], "costType"));
	float smoothFactor = nrhs > 8 ? getCastScalarChecked<float>(prhs[8], "smoothFactor") : 0.0f;
	int smoothIter = nrhs > 9 ? getCastScalarChecked<int>(prhs[9], "smoothIter") : 0;

	bool scalarRadii = initRs.nx == 1 || initRs.ny == 1;
	int nr = scalarRadii ? initRs.numElem() : initRs.nz;
	int nmesh = centers.nx;
	ensureOrError(nr == 1 || nr == nmesh, "Must provide a single initR or one for every center");
	int nthreads = std::min(getMaxCompThreads(), nmesh);
	int meshesPerThread = (nmesh + (nmesh % nthreads == 0 ? 0 : nthreads)) / nthreads;

	std::vector<float *> vertData;
	std::vector<int *> faceData;

	// Make worker function to process subset
	auto workerFunc = [&](int begin, int end) {
		for (int i = begin; i < end; ++i) {
			float initR = scalarRadii ? initRs[i] : 1.0f;
			Vec3f center(centers.at(i, 0, 0), centers.at(i, 1, 0), centers.at(i, 2, 0));
			ManifoldMesh init = SubdividedIcosahedron(center, initR, initSubDiv);
			if (!scalarRadii) {
				Mat3x3f A(
					Vec3f(initRs.at(0, 0, i), initRs.at(1, 0, i), initRs.at(2, 0, i)),
					Vec3f(initRs.at(0, 1, i), initRs.at(1, 1, i), initRs.at(2, 1, i)),
					Vec3f(initRs.at(0, 2, i), initRs.at(1, 2, i), initRs.at(2, 2, i))
				);
				for (auto& v : init.vertices) {
					v.pos = center + A*(v.pos - center);
				}
			}
			ManifoldMesh mesh = surfaceCut(cost, std::move(init), numSamples, sampleStep, maxDiff, costType);

			mesh.taubinSmooth(smoothFactor, 1.04f * smoothFactor, smoothIter);

			const size_t numVerts = mesh.vertices.size();
			const size_t numFaces = mesh.faces.size();

			for (const auto& v : mesh.vertices) {
				vertData[i][v.self + 0 * numVerts] = v.pos[0];
				vertData[i][v.self + 1 * numVerts] = v.pos[1];
				vertData[i][v.self + 2 * numVerts] = v.pos[2];
			}

			for (const auto& f : mesh.faces) {
				// Add one since MATLAB uses 1-indexing
				faceData[i][f.self + 0 * numFaces] = mesh.edges[f.edge].vert + 1;
				faceData[i][f.self + 1 * numFaces] = mesh.edges[mesh.next(f.edge)].vert + 1;
				faceData[i][f.self + 2 * numFaces] = mesh.edges[mesh.next(mesh.next(f.edge))].vert + 1;
			}
		}
	};

	// Allocate all needed memory
	mxArray *vertCell = mxCreateCellMatrix(1, nmesh);
	mxArray *faceCell = mxCreateCellMatrix(1, nmesh);
	vertData.reserve(nmesh);
	faceData.reserve(nmesh);
	const size_t triNum = std::pow(4, initSubDiv);
	const size_t numFaces = 20 * triNum;
	const size_t numVerts = 10 * triNum + 2;
	for (int i = 0; i < nmesh; ++i) {
		mxArray *vertices = mxCreateNumericMatrix(numVerts, 3, mxSINGLE_CLASS, mxREAL);
		mxArray *faces = mxCreateNumericMatrix(numFaces, 3, mxINT32_CLASS, mxREAL);
		vertData.push_back(static_cast<float *>(mxGetData(vertices)));
		faceData.push_back(static_cast<int *>(mxGetData(faces)));

		mxSetCell(vertCell, i, vertices);
		mxSetCell(faceCell, i, faces);
	}

	// Compute all cuts
	if (nthreads == 1) {
		// Just run on this thread
		workerFunc(0, nmesh);
	} else {
		std::vector<std::thread> threads;
		// Split processing among threads
		for (int i = 0; i < nthreads; ++i) {
			int begin = i * meshesPerThread;
			int end = std::min((i + 1) * meshesPerThread, nmesh);
			threads.push_back(std::thread(workerFunc, begin, end));
		}
		// Wait for all threads to finish
		for (auto& th : threads) {
			if (th.joinable()) {
				th.join();
			}
		}
	}

	plhs[0] = faceCell;
	plhs[1] = vertCell;
}
