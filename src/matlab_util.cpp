#include "matlab_util.h"
#include <algorithm>
#include <sstream>

bool isSize(const mxArray *a, std::initializer_list<int> size)
{
	size_t dimA = mxGetNumberOfDimensions(a);
	int i = 0;
	for (int s : size) {
		if (s != -1 && (i >= dimA || s != mxGetDimensions(a)[i])) {
			return false;
		}
		++i;
	}
	return true;
}

void ensureSize(const mxArray* a, std::initializer_list<int> size, const std::string& errPrefixStr, 
    const std::string& wildcard)
{
    if (!isSize(a, size)) {
        const int numWildcards = std::count(size.begin(), size.end(), -1);
        int crntWildcard = 1;
        // Build string for expected size
        // This is likely not the fastest way to construct this, but we are aborting right after
        std::ostringstream oss;
        int i = 0;
        for (int s : size) {
            if (s != -1) {
                oss << s;
            } else {
                oss << wildcard;
                if (numWildcards > 1) {
                    oss << crntWildcard;
                    ++crntWildcard;
                }
            }
            if (i < size.size() - 1) {
                oss << " x ";
            }
            ++i;
        }
        std::string expectedSizeStr = oss.str();
        
        // Build string for actual size
        oss.str(""); // Clear stringstream
        const int numDim = mxGetNumberOfDimensions(a);
        for (i = 0; i < numDim; ++i) {
            oss << mxGetDimensions(a)[i];
            if (i < numDim - 1) {
                oss << " x ";
            }
        }
        std::string actualSizeStr = oss.str();
        abortWithMsg("%s must be an %s array (was: %s)", errPrefixStr, expectedSizeStr, actualSizeStr);
    }
}

void ensureArgCount(int narg, int num)
{
    ensureOrError(narg == num, "Must supply %d inputs", num);
}

void ensureArgRange(int narg, int min, int max)
{
    if (min == max) {
        ensureArgCount(narg, min);
    } else {
        ensureOrError(min <= narg && narg <= max, "Must supply between %d and %d inputs", min, max);
    }
}

int getMaxCompThreads()
{
	mxArray *matlabCallOut[1] = { 0 };
	mxArray *matlabCallIn[1] = { 0 };
	mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
	double *Nthreadsd = mxGetPr(matlabCallOut[0]);
	return (int)Nthreadsd[0];
}
