#ifndef VOLUME_H__
#define VOLUME_H__

#include <cstdint>
#include <memory>

#include <GEL/CGLA/Vec3i.h>
#include <GEL/CGLA/Vec3f.h>

using namespace CGLA;

/** Simple structure to hold volumetric data of type Ty */
template <class Ty>
struct Volume {
    typedef Ty DataType;

    std::shared_ptr<Ty> data;
    size_t nx;
    size_t ny;
    size_t nz;

    explicit Volume() :
        data(),
        nx(0),
        ny(0),
        nz(0) {}
    explicit Volume(size_t nx, size_t ny, size_t nz) noexcept;
    Volume(const Volume&) noexcept;
    Volume(Volume&& other) noexcept;

    Volume& operator=(const Volume&) noexcept;
    Volume& operator=(Volume&& other) noexcept;

    void alloc();
    void alloc(size_t nelem);
    void alloc(size_t nx, size_t ny, size_t nz);

    void clear();

	void copy(const Volume<Ty>& src);
	void copy(Volume<Ty>&& src) noexcept;

    size_t numElem() const noexcept;

    Ty interp(float x, float y, float z) const;
    Ty interp(const Vec3f& p) const;

    size_t idx(size_t x, size_t y, size_t z) const noexcept;
    Vec3i pos(size_t idx) const noexcept;
	bool contains(const Vec3f& p) const noexcept;

    Ty& at(size_t i);
    const Ty& at(size_t i) const;
    Ty& at(size_t x, size_t y, size_t z);
    const Ty& at(size_t x, size_t y, size_t z) const;
    Ty& at(const Vec3i& p);
    const Ty& at(const Vec3i& p) const;

    Ty& operator[](size_t i);
    const Ty& operator[](size_t i) const;
    Ty& operator[](const Vec3i& p);
    const Ty& operator[](const Vec3i& p) const;

    void print() const;
    void saveToBin(const std::string& fname) const;
    void loadFromBin(const std::string& fname);
};

template <class Ty>
constexpr uint32_t getTypeCode();

#include <volume.inl>

#endif // VOLUME_H__
