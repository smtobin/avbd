#pragma once

#include "common/mesh/ParticleMesh.hpp"

class ParticleTetMesh : public ParticleMesh
{
protected:
    /** Tetrahedral elements - 4-vectors of indices in the vertices list
     * (they do not correspond to indices in the particle pool - must go through the vertices list)
     */
    TombstoneVector<Vec4u> _elements;

public:
    ParticleTetMesh() = default;
    
    /** Constructs a tetrahedral mesh from a set of vertices, faces, and elements.
     * This is usually done using the helper methods in the MeshUtils library.
     */
    ParticleTetMesh(ParticlePool& pool, const std::vector<Vec3r>& vertices, const std::vector<Vec3u>& faces, const std::vector<Vec4u>& elements);

    /** Returns a const-reference to the elements of the mesh. */
    const TombstoneVector<Vec4u>& elements() const { return _elements; }

    /** Returns the number of elements in the mesh. */
    int numElements() const { return _elements.size(); }

    /** Returns a single element as an Eigen 4-vector, given the element index. */
    const Vec4u& element(int index) const { return _elements.at(index); }

    /** Returns whether or not the index corresponds to a valid element. */
    bool elementValid(int index) const { return _elements.indexValid(index); }

    /** Returns the current volume of the specified element. */
    Real elementVolume(int index) const;

    /** Returns the current centroid of the element. */
    Vec3r elementCentroid(int elem_index) const;
};