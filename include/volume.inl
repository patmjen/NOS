#include "volume.h"
#include <fstream>
#include <cassert>
#include <iostream>
#include <type_traits>

template <class Ty>
Volume<Ty>::Volume(size_t nx, size_t ny, size_t nz) noexcept :
    data(),
    nx(nx),
    ny(ny),
    nz(nz)
{}

template <class Ty>
Volume<Ty>::Volume(const Volume& other) noexcept :
    data(other.data),
    nx(other.nx),
    ny(other.ny),
    nz(other.nz)
{}

template <class Ty>
Volume<Ty>::Volume(Volume&& other) noexcept :
    data(std::move(other.data)),
    nx(other.nx),
    ny(other.ny),
    nz(other.nz)
{
	other.nx = 0;
	other.ny = 0;
	other.nz = 0;
}

template <class Ty>
Volume<Ty>& Volume<Ty>::operator=(const Volume& other) noexcept
{
    if (this != &other) {
        data = other.data;
        nx = other.nx;
        nz = other.nz;
        ny = other.ny;
    }
    return *this;
}

template <class Ty>
Volume<Ty>& Volume<Ty>::operator=(Volume&& other) noexcept
{
    if (this != &other) {
        data = std::move(other.data);
        nx = other.nx;
        ny = other.ny;
        nz = other.nz;

		other.nx = 0;
		other.ny = 0;
		other.nz = 0;
    }
    return *this;
}

template <class Ty>
void Volume<Ty>::alloc()
{
    alloc(nx, ny, nz);
}

template <class Ty>
void Volume<Ty>::alloc(size_t nelem)
{
    // Need to supply custom deleter since shared_ptr does not support
    // arrays prior to C++17
    data = std::shared_ptr<Ty>(new Ty[nelem], std::default_delete<Ty[]>());
}

template <class Ty>
void Volume<Ty>::alloc(size_t nx, size_t ny, size_t nz)
{
    alloc(nx * ny * nz);
}

template <class Ty>
void Volume<Ty>::clear() {
    nx = 0;
    ny = 0;
    nz = 0;
    data = nullptr;
}

template <class Ty>
void Volume<Ty>::copy(const Volume<Ty>& src)
{
	if (this != &src) {
		nx = src.nx;
		ny = src.ny;
		nz = src.nz;
		alloc();
		std::memcpy(data.get(), src.get(), nx * ny * nz * sizeof(Ty));
	}
}

template <class Ty>
void Volume<Ty>::copy(Volume<Ty>&& src) noexcept
{
	if (this != &src) {
		// If src is an r-value we can just move from it instead of doing a deep copy
		*this = std::move(src);
	}
}

template <class Ty>
size_t Volume<Ty>::numElem() const noexcept
{
    return nx * ny * nz;
}

template <class Ty>
Ty Volume<Ty>::interp(float x, float y, float z) const
{
    float fx = x - floorf(x);
    float fy = y - floorf(y);
    float fz = z - floorf(z);
    size_t x0 = static_cast<size_t>(floorf(x));
    size_t y0 = static_cast<size_t>(floorf(y));
    size_t z0 = static_cast<size_t>(floorf(z));
    // Handle special case where query point is on boundary by just re-using first point
    size_t x1 = (x0 == nx - 1 && x <= nx - 1) ? x0 : x0 + 1;
    size_t y1 = (y0 == ny - 1 && y <= ny - 1) ? y0 : y0 + 1;
    size_t z1 = (z0 == nz - 1 && z <= nz - 1) ? z0 : z0 + 1;

    // Do interpolation
    Ty c000 = at(x0, y0, z0);
    Ty c001 = at(x0, y0, z1);
    Ty c010 = at(x0, y1, z0);
    Ty c011 = at(x0, y1, z1);
    Ty c100 = at(x1, y0, z0);
    Ty c101 = at(x1, y0, z1);
    Ty c110 = at(x1, y1, z0);
    Ty c111 = at(x1, y1, z1);
    Ty c00 = (1.0f - fx) * c000 + fx * c100;
    Ty c01 = (1.0f - fx) * c001 + fx * c101;
    Ty c10 = (1.0f - fx) * c010 + fx * c110;
    Ty c11 = (1.0f - fx) * c011 + fx * c111;
    Ty c0 = (1.0f - fy) * c00 + fy * c10;
    Ty c1 = (1.0f - fy) * c01 + fy * c11;
    return (1.0f - fz) * c0 + fz * c1;
}

template <class Ty>
Ty Volume<Ty>::interp(const Vec3f& p) const
{
    return interp(p[0], p[1], p[2]);
}

template <class Ty>
size_t Volume<Ty>::idx(size_t x, size_t y, size_t z) const noexcept
{
    return x + y * nx + z * nx * ny;
}

template <class Ty>
Vec3i Volume<Ty>::pos(size_t idx) const noexcept
{
    return Vec3i(
        idx % nx,
        (idx / nx) % ny,
        idx / (nx * ny)
    );
}

template <class Ty>
bool Volume<Ty>::contains(const Vec3f& p) const noexcept
{
	return 0 <= p[0] && 0 <= p[1] && 0 <= p[2] &&
		p[0] <= nx - 1 && p[1] <= ny - 1 && p[2] <= nz - 1;
}

template <class Ty>
Ty& Volume<Ty>::at(size_t i)
{
    assert(i < nx * ny * nz);
    return data.get()[i];
}

template <class Ty>
const Ty& Volume<Ty>::at(size_t i) const
{
    assert(i < nx * ny * nz);
    return data.get()[i];
}

template <class Ty>
Ty& Volume<Ty>::at(size_t x, size_t y, size_t z)
{
    return at(idx(x, y, z));
}

template <class Ty>
const Ty& Volume<Ty>::at(size_t x, size_t y, size_t z) const
{
    return at(idx(x, y, z));
}

template <class Ty>
Ty& Volume<Ty>::at(const Vec3i& p)
{
    assert(p[0] >= 0 && p[1] >= 0 && p[2] >= 0);
    return at(p[0], p[1], p[2]);
}

template <class Ty>
const Ty& Volume<Ty>::at(const Vec3i& p) const
{
    assert(p[0] >= 0 && p[1] >= 0 && p[2] >= 0);
    return at(p[0], p[1], p[2]);
}

template <class Ty>
Ty& Volume<Ty>::operator[](size_t i)
{
    return at(i);
}

template <class Ty>
const Ty& Volume<Ty>::operator[](size_t i) const
{
    return at(i);
}

template <class Ty>
Ty& Volume<Ty>::operator[](const Vec3i &p)
{
    return at(p);
}

template <class Ty>
const Ty& Volume<Ty>::operator[](const Vec3i& p) const
{
    return at(p);
}

template <class Ty>
void Volume<Ty>::print() const
{
    // Print the volume like NumPy does in Python
    std::cout << '[';
    for (int z = 0; z < nz; ++z) {
        if (z > 0) {
            std::cout << ' ';
        }
        std::cout << '[';
        for (int y = 0; y < ny; ++y) {
            std::cout << "[ ";
            for (int x = 0; x < nx; ++x) {
                std::cout << at(x, y, z);
                if (x + 1 < nx) {
                    std::cout << ' ';
                }
            }
            std::cout << " ]";
            if (y + 1 < ny) {
                std::cout << "\n  ";
            }
        }
        std::cout << ']';
        if (z + 1 < nz) {
            std::cout << "\n\n";
        }
    }
    std::cout << "]\n";
}

template <class Ty>
void Volume<Ty>::saveToBin(const std::string& fname) const
{
    // Write volume to simple binary file format with format
    // 1. (uint32_t) Type, see getTypeCode
    // 2. (uint32_t) nx
    // 3. (uint32_t) ny
    // 4. (uint32_t) nz
    // 5. (Ty) data[0] (Ty) data[1] (Ty) data[2]...
    uint32_t typeU32 = getTypeCode<Ty>();
    uint32_t nxU32 = nx;
    uint32_t nyU32 = ny;
    uint32_t nzU32 = nz;

    std::fstream file(fname, std::ios::out | std::ios::binary);
    assert(file); // Ensure file was opened

    // Write header
    file.write((char *)&typeU32, sizeof(uint32_t));
    file.write((char *)&nxU32, sizeof(uint32_t));
    file.write((char *)&nyU32, sizeof(uint32_t));
    file.write((char *)&nzU32, sizeof(uint32_t));

    // Write data
    file.write((char *)data.get(), numElem() * sizeof(Ty));
}

template <class Ty>
void Volume<Ty>::loadFromBin(const std::string& fname)
{
    uint32_t typeCode;
    uint32_t nxU32;
    uint32_t nyU32;
    uint32_t nzU32;

    std::fstream file(fname, std::ios::in | std::ios::binary);
    assert(file); // Ensure file was opened

    // Read type code and ensure it is the same as this volume
    file.read((char *)&typeCode, sizeof(uint32_t));
    assert(typeCode == getTypeCode<Ty>()); // Must load data of right type

    clear(); // Prepare for new data

    // Read new size and alloc space
    file.read((char *)&nxU32, sizeof(uint32_t));
    file.read((char *)&nyU32, sizeof(uint32_t));
    file.read((char *)&nzU32, sizeof(uint32_t));
    nx = nxU32;
    ny = nyU32;
    nz = nzU32;
    alloc();

    // Read data
    file.read((char *)data.get(), numElem() * sizeof(Ty));
}

template <class Ty>
constexpr uint32_t getTypeCode()
{
	// Codes are from Table 1-1 in MAT - File Format
	// URL: https://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf

    // We cannot use switch so have to revert to this
    if (std::is_same<Ty, int8_t>::value) {
        return 1;
    } else if (std::is_same<Ty, uint8_t>::value) {
        return 2;
    } else if (std::is_same<Ty, int16_t>::value) {
        return 3;
    } else if (std::is_same<Ty, uint16_t>::value) {
        return 4;
    } else if (std::is_same<Ty, int32_t>::value) {
        return 5;
    } else if (std::is_same<Ty, uint32_t>::value) {
        return 6;
    } else if (std::is_same<Ty, float>::value) {
        return 7;
    } else if (std::is_same<Ty, double>::value) {
        return 9;
    } else if (std::is_same<Ty, int64_t>::value) {
        return 12;
    } else if (std::is_same<Ty, uint64_t>::value) {
        return 13;
    }
}