#pragma once

#include "common/common.hpp"
#include "common/ParticlePool.hpp"

#include "common/TombstoneVector.hpp"

class ParticleMesh
{
protected:
    /** Pointer to the particle memory pool */
    ParticlePool* _particle_pool;

    /** Vertices - indices in the particle pool.
     */
    TombstoneVector<unsigned> _vertices;

    /** Faces - 3-vectors of indices in the vertices list
     * (they do not correspond to indices in the particle pool - must go through the vertices list)
     */
    TombstoneVector<Vec3u> _faces;

    /** The size of the mesh in its unrotated state - i.e. the size of the oriented bounding box */
    // Vec3r _unrotated_size_xyz;

    /** The current location of the original (0,0,0) point in mesh space. */
    Vec3r _mesh_origin;

public:
    ParticleMesh() = default;

    /** Constructs a mesh from a set of vertices and faces.
     * This is usually done using helper methods in the MeshUtils library.
     */
    ParticleMesh(ParticlePool& pool, const std::vector<Vec3r> &vertices, const std::vector<Vec3u> &faces);

    /** Returns a const-reference to the vertices of the mesh. */
    const TombstoneVector<unsigned>& vertices() const { return _vertices; }
    /** Returns a const-reference to the faces of the mesh. */
    const TombstoneVector<Vec3u>& faces() const { return _faces; }

    /** Number of verticees in the mesh. */
    int numVertices() const { return _vertices.size(); }
    /** Number of faces in the mesh. */
    int numFaces() const { return _faces.size(); }

    /** Essentially "sets up" the mesh - treats the current state as the initial, undeformed state of the mesh.
     * This should be called after performing the initial translations and rotations setting up the mesh.
     */
    virtual void setCurrentStateAsUndeformedState();

    /** Returns a single vertex as an Eigen 3-vector, given the vertex index.
     * This assumes that the index used is a valid index (i.e. the vertex we are trying to access has not been removed).
     */
    const Vec3r& vertex(int index) const { return _particle_pool->positions[_vertices[index]]; }
    const Vec3r& previousVertex(int index) const { return _particle_pool->previous_positions[_vertices[index]]; }

    /** Returns whether not the index corresponds to a valid vertex. */
    bool vertexValid(int index) const { return _particle_pool->active[_vertices[index]]; }

    /** Sets the vertex at the specified to a new position. */
    void setVertex(int index, const Vec3r &new_pos) { _particle_pool->positions[_vertices.at(index)] = new_pos; }

    /** Returns the surface face at the given index.
     * This assumes that the index used is a valid index (i.e. the face we are trying to access has not been removed).
     */
    const Vec3u& face(int index) const { return _faces.at(index); }

    /** Returns whether or not the index corresponds to a valid face. */
    bool faceValid(int index) const { return _faces.indexValid(index); }

    /** Returns the coordinates of the mesh origin.
     * i.e. where the "origin" of the mesh (the (0,0,0) point when the mesh is first loaded) is currently at
     */
    const Vec3r& meshOrigin() const { return _mesh_origin; }

    /** Returns the unrotated size of the mesh.
     * This does not change when the mesh rotates - only when the mesh is resized.
     */
    // Vec3r unrotatedSize() const { return _unrotated_size_xyz; }

    /** Scales the mesh according to the given scaling values in the x, y, and z directions. */
    void scale(const Vec3r& scaling);

    /** Moves each vertex in the mesh by the same amount. */
    void moveTogether(const Vec3r &delta);

    /** Rotates the mesh according to a vector of (x angle, y angle, and z angle) Euler angles around some point p.
     * Usees the XYZ Euler angle convention - rotates x degrees about x-axis, then y degrees about y-axis, and then z degrees about z-axis
     * @param p : the point around which to rotate the mesh
     * @param xyz_angles : a 3-vector corresponding to successive rotation angles
     */
    // void rotateAbout(const Vec3r &p, const Vec3r &xyz_angles);

    /** Rotates the mesh using a 3x3 rotation matrix around some point p.
     * @param p : the point around which to rotate the mesh
     * @param rot_mat : the rotation matrix used to rotate the vertices
     */
    void rotateAbout(const Vec3r &p, const Mat3r &rot_mat);

    /** Computes the current total mass, center of mass, and moment of inertia tensor (about center of mass) for the mesh, given a density.
     * Uses the algorithm described here: http://number-none.com/blow/inertia/index.html
     * @param density : the density to be used in the calculation (DEFAULT = 1). If omitted, the "mass" returned is actually the volume of the mesh.
     * @returns these quantities as a 3-tuple: (mass, center-of-mass, moment-of-inertia)
     */
    std::tuple<Real, Vec3r, Mat3r> massProperties(Real density = 1.0) const;

    /** Computes the current center of mass for the mesh.
     * Uses the same algorithm as massProperties(), but without calculating the moment of inertia.
     */
    Vec3r massCenter() const;

    /** Writes the mesh to .obj file.
     * Only writes vertices and faces.
     */
    void writeMeshToObjFile(const std::string &filename) const;

    /** Checks if a point p is inside the mesh. Uses winding number approach. (O(n) computation) */
    bool isInside(const Vec3r& p) const;
};