#ifndef MATLAB_UTIL_H__
#define MATLAB_UTIL_H__

#include <string>
#include <memory>
#include <typeinfo>
#include <cstdint>
#include <initializer_list>
#include <type_traits>

#include "mex.h"
#include "matrix.h"

#include <GEL/CGLA/ArithVec.h>

#include "volume.h"

#define abortWithMsg(...) do { \
    mexErrMsgIdAndTxt("surfseg:internal", __VA_ARGS__); \
} while(false)

#define ensureOrError(expr, ...) do { \
	if (!(expr)) { \
		abortWithMsg(__VA_ARGS__); \
	} \
} while(false)

#ifdef __cplusplus
extern "C" bool utIsInterruptPending();
#else
extern bool utIsInterruptPending();
#endif

template <class Ty>
inline mxClassID typeClassId()
{
	// We cannot use switch so have to revert to this
	if (std::is_same<Ty, mxLogical>::value) {
		return mxLOGICAL_CLASS;
	} else if (std::is_same<Ty, mxInt8>::value) {
		return mxINT8_CLASS;
	} else if (std::is_same<Ty, mxUint8>::value) {
		return mxUINT8_CLASS;
	} else if (std::is_same<Ty, mxInt16>::value) {
		return mxINT16_CLASS;
	} else if (std::is_same<Ty, mxUint16>::value) {
		return mxUINT16_CLASS;
	} else if (std::is_same<Ty, mxInt32>::value) {
		return mxINT32_CLASS;
	} else if (std::is_same<Ty, mxUint32>::value) {
		return mxUINT32_CLASS;
	} else if (std::is_same<Ty, mxSingle>::value) {
		return mxSINGLE_CLASS;
	} else if (std::is_same<Ty, mxDouble>::value) {
		return mxDOUBLE_CLASS;
	} else if (std::is_same<Ty, mxInt64>::value) {
		return mxINT64_CLASS;
	} else if (std::is_same<Ty, mxUint64>::value) {
		return mxUINT64_CLASS;
	} else {
		return mxUNKNOWN_CLASS;
	}
}

inline const char *getClassNameFromId(mxClassID id)
{
	switch (id) {
	case mxLOGICAL_CLASS:
		return "logical";
	case mxINT8_CLASS:
		return "int8";
	case mxUINT8_CLASS:
		return "uint8";
	case mxINT16_CLASS:
		return "int16";
	case mxUINT16_CLASS:
		return "uint16";
	case mxINT32_CLASS:
		return "int32";
	case mxUINT32_CLASS:
		return "uint32";
	case mxSINGLE_CLASS:
		return "single";
	case mxDOUBLE_CLASS:
		return "double";
	case mxINT64_CLASS:
		return "int64";
	case mxUINT64_CLASS:
		return "uint64";
	default:
		return "unkown";
	}
}

template <class Ty>
inline Ty derefAndCast(const void *ptr, mxClassID id)
{
	switch (id) {
	case mxLOGICAL_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxLogical *>(ptr));
	case mxINT8_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxInt8 *>(ptr));
	case mxUINT8_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxUint8 *>(ptr));
	case mxINT16_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxInt16 *>(ptr));
	case mxUINT16_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxUint16 *>(ptr));
	case mxINT32_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxInt32 *>(ptr));
	case mxUINT32_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxUint32 *>(ptr));
	case mxSINGLE_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxSingle *>(ptr));
	case mxDOUBLE_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxDouble *>(ptr));
	case mxINT64_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxInt64 *>(ptr));
	case mxUINT64_CLASS:
		return static_cast<Ty>(*reinterpret_cast<const mxUint64 *>(ptr));
	default:
		abortWithMsg("derefAndCast: Unknown class ID");
		return Ty(); // NOTE: This will never run as mexErrMsgIdAndTxt aborts execution
	}
}

inline bool isVector(const mxArray *a)
{
	return mxGetNumberOfDimensions(a) == 2 && (mxGetN(a) == 1 || mxGetM(a) == 1);
}

inline bool isMatrix(const mxArray *a)
{
	return mxGetNumberOfDimensions(a) == 2;
}

bool isSize(const mxArray *a, std::initializer_list<int> size);

void ensureSize(const mxArray* a, std::initializer_list<int> size,
    const std::string& errPrefixStr = "Value", const std::string& wildcard = "N");

template <class Ty>
inline void ensureMatchingClass(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	const char *errPrefix = errPrefixStr.c_str();
	ensureOrError(mxGetClassID(a) == typeClassId<Ty>(), "%s must be of class %s (was: %s)",
		errPrefix, getClassNameFromId(typeClassId<Ty>()), mxGetClassName(a));
}

void ensureArgCount(int narg, int num);

void ensureArgRange(int narg, int min, int max);

template <class Ty>
inline Ty getCastScalarChecked(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	const char *errPrefix = errPrefixStr.c_str();

	ensureOrError(!mxIsComplex(a), "%s must be real", errPrefix);
	ensureOrError(mxIsScalar(a), "%s must be a scalar", errPrefix);

	return static_cast<Ty>(mxGetScalar(a));
}

template <class Ty, unsigned int N>
inline std::array<Ty, N> getVectorChecked(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	mxClassID classId = mxGetClassID(a);
	const char *errPrefix = errPrefixStr.c_str();

	ensureOrError(!mxIsComplex(a), "%s must be real", errPrefix);
	ensureMatchingClass<Ty>(a, errPrefixStr);
	ensureOrError(isVector(a), "%s must be a vector", errPrefix);
	ensureOrError(mxGetNumberOfElements(a) == N, "%s must have length %d", errPrefix, N);

	std::array<Ty, N> out;
	const Ty *data = static_cast<Ty *>(mxGetData(a));
	for (int i = 0; i < N; ++i) {
		out[i] = data[i];
	}
	return out;
}

template <class VecTy>
inline VecTy getVectorChecked(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	typedef typename VecTy::ScalarType Ty;
	const size_t N = VecTy::get_dim();
	mxClassID classId = mxGetClassID(a);
	const char *errPrefix = errPrefixStr.c_str();

	ensureOrError(!mxIsComplex(a), "%s must be real", errPrefix);
	ensureMatchingClass<Ty>(a, errPrefixStr);
	ensureOrError(isVector(a), "%s must be a vector", errPrefix);
	ensureOrError(mxGetNumberOfElements(a) == N, "%s must have length %d", errPrefix, N);

	VecTy out;
	const Ty *data = static_cast<Ty *>(mxGetData(a));
	for (int i = 0; i < N; ++i) {
		out[i] = data[i];
	}
	return out;
}

template <class VecTy>
inline VecTy getCastVectorChecked(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	typedef typename VecTy::ScalarType Ty;
	const size_t N = VecTy::get_dim();
	mxClassID classId = mxGetClassID(a);
	const char *errPrefix = errPrefixStr.c_str();

	ensureOrError(mxIsNumeric(a), "%s must be numeric", errPrefix);
	ensureOrError(!mxIsComplex(a), "%s must be real", errPrefix);
	ensureOrError(isVector(a), "%s must be a vector", errPrefix);
	ensureOrError(mxGetNumberOfElements(a) == N, "%s must have length %d", errPrefix, N);

	VecTy out;
	const std::uint8_t *data = static_cast<const std::uint8_t *>(mxGetData(a));
	for (int i = 0; i < N; ++i) {
		out[i] = derefAndCast<Ty>(data, classId);
		data += mxGetElementSize(a);
	}
	return out;
}

template <class Ty>
inline std::shared_ptr<Ty> getCastSharedPtrChecked(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	mxClassID classId = mxGetClassID(a);
	const char *errPrefix = errPrefixStr.c_str();

	ensureOrError(mxIsNumeric(a), "%s must be numeric", errPrefix);
	ensureOrError(!mxIsComplex(a), "%s must be real", errPrefix);

	if (classId == typeClassId<Ty>()) {
		// Type of MATLAB array already matches type, so just use it
		return std::shared_ptr<Ty>(static_cast<Ty *>(mxGetData(a)), [](auto) {});
	} else {
		// Type does *not* match so alloc an array to hold and cast all values
		size_t n = mxGetNumberOfElements(a);
		std::shared_ptr<Ty> out(static_cast<Ty *>(mxMalloc(n * sizeof(Ty))), [](auto ptr) { mxFree(ptr); });
		Ty *outPtr = out.get();
		// Need to use uint8_t pointer so we can do pointer arithmetic
		const std::uint8_t *data = static_cast<const std::uint8_t *>(mxGetData(a));
		for (int i = 0; i < n; ++i) {
			outPtr[i] = derefAndCast<Ty>(data, classId);
			data += mxGetElementSize(a);
		}
		return out;
	}
}

template <class Ty>
inline Volume<Ty> getVolumeChecked(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	mxClassID classId = mxGetClassID(a);
	const char *errPrefix = errPrefixStr.c_str();

	ensureOrError(!mxIsComplex(a), "%s must be real", errPrefix);
	ensureMatchingClass<Ty>(a, errPrefixStr);
	ensureOrError(mxGetNumberOfDimensions(a) <= 3, "%s have must at most 3 dimensions", errPrefix);

	const mwSize *dims = mxGetDimensions(a);
	Volume<Ty> out;
	out.nx = dims[0];
	out.ny = dims[1];
	out.nz = mxGetNumberOfDimensions(a) > 2 ? dims[2] : 1;
	// Make a shared pointer with no deleter, since the memory comes from MATLAB
	out.data = std::shared_ptr<Ty>(static_cast<Ty *>(mxGetData(a)), [](auto) {});
	return out;
}

template <class Ty>
inline Volume<Ty> getCastVolumeChecked(const mxArray *a, const std::string& errPrefixStr = "Value")
{
	mxClassID classId = mxGetClassID(a);
	const char *errPrefix = errPrefixStr.c_str();

	ensureOrError(!mxIsComplex(a), "%s must be real", errPrefix);
	ensureOrError(mxGetNumberOfDimensions(a) <= 3, "%s have must at most 3 dimensions", errPrefix);

	const mwSize *dims = mxGetDimensions(a);
	Volume<Ty> out;
	out.nx = dims[0];
	out.ny = dims[1];
	out.nz = mxGetNumberOfDimensions(a) > 2 ? dims[2] : 1;
	out.data = getCastSharedPtrChecked<Ty>(a, errPrefixStr);
	return out;
}

int getMaxCompThreads();

#endif // MATLAB_UTIL_H__
