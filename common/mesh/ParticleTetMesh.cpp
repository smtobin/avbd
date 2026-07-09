#include "common/mesh/ParticleTetMesh.hpp"

ParticleTetMesh::ParticleTetMesh(ParticlePool& pool, const std::vector<Vec3r>& vertices, const std::vector<Vec3u>& faces, const std::vector<Vec4u>& elements)
    : ParticleMesh(pool, vertices, faces), _elements(elements)
{

}

Real ParticleTetMesh::elementVolume(int index) const
{
    const Vec4u& elem = element(index);
    const Vec3r& v1 = vertex(elem[0]);
    const Vec3r& v2 = vertex(elem[1]);
    const Vec3r& v3 = vertex(elem[2]);
    const Vec3r& v4 = vertex(elem[3]);

    Mat3r X;
    X.col(0) = (v1 - v4);
    X.col(1) = (v2 - v4);
    X.col(2) = (v3 - v4);

    return std::abs(X.determinant() / 6.0);
}

Vec3r ParticleTetMesh::elementCentroid(int element_index) const
{
    const Vec4u& elem_verts = element(element_index);
    Vec3r sum = vertex(elem_verts[0]) + vertex(elem_verts[1]) + vertex(elem_verts[2]) + vertex(elem_verts[3]);
    return sum/4.0;
}