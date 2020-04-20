#ifndef SUBDIVIDED_ICOSAHEDRON_H__
#define SUBDIVIDED_ICOSAHEDRON_H__

#include <GEL/CGLA/Vec3f.h>
#include <GEL/CGLA/Mat3x3f.h>

#include "manifold_mesh.h"

using namespace CGLA;

class SubdividedIcosahedron : public ManifoldMesh {
    Vec3f center_;
    float r_;
    int divisionLevel_;

public:
    explicit SubdividedIcosahedron() :
        ManifoldMesh(),
        center_(0.0f),
        r_(0.0f),
        divisionLevel_(1) {};
    explicit SubdividedIcosahedron(float r);
    explicit SubdividedIcosahedron(const Vec3f& center, float r, int divLvl = 1);
    explicit SubdividedIcosahedron(const Mat3x3f& R);
    explicit SubdividedIcosahedron(const Vec3f& center, const Mat3x3f& R, int divLvl = 1);

    void buildIcosahedron();
    void buildIcosahedron(const Vec3f& center, float r);

    void rescale(float newR) noexcept;
    void move(const Vec3f& newCenter) noexcept;
    void moveAndRescale(const Vec3f& newCenter, float newR) noexcept;

    float r() const noexcept;
    const Vec3f& center() const noexcept;

    void subdivide(int numDivide = 1);

private:
    void singleSubdivide_();
};

#endif // SUBDIVIDED_ICOSAHEDRON_H__
