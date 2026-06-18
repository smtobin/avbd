#include "common/mesh/ParticleMesh.hpp"

#include <fstream>

ParticleMesh::ParticleMesh(ParticlePool& pool, const std::vector<Vec3r>& vertices, const std::vector<Vec3u>& faces)
    : _particle_pool(&pool), _faces(faces)
{
    // create particles from vertices list
    for (const auto& vert : vertices)
    {
        unsigned idx = _particle_pool->addParticle(vert, 0);
        _vertices.push_back(idx);
    }

    _mesh_origin = Vec3r::Zero();
}

void ParticleMesh::scale(const Vec3r& scaling)
{
    for (const auto& v : _vertices)
    {
        _particle_pool->positions[v][0] *= scaling[0];
        _particle_pool->positions[v][1] *= scaling[1];
        _particle_pool->positions[v][2] *= scaling[2];
    }

    // for (auto& v : _initial_vertices)
    // {
    //     v[0] *= scaling[0];
    //     v[1] *= scaling[1];
    //     v[2] *= scaling[2];
    // }

    _mesh_origin[0] *= scaling[0];
    _mesh_origin[1] *= scaling[1];
    _mesh_origin[2] *= scaling[2];

    // scale the unrotated size
    // _unrotated_size_xyz[0] *= scaling[0];
    // _unrotated_size_xyz[1] *= scaling[1];
    // _unrotated_size_xyz[2] *= scaling[2];
}

void ParticleMesh::moveTogether(const Vec3r& delta)
{
    for (const auto& v : _vertices)
        _particle_pool->positions[v] += delta;

    // for (auto& v : _initial_vertices)
    //     v += delta;
        
    _mesh_origin += delta;
}

void ParticleMesh::rotateAbout(const Vec3r& p, const Mat3r& rot_mat)
{
    moveTogether(-p);
    for (auto& v : _vertices)
        _particle_pool->positions[v] = rot_mat * _particle_pool->positions[v];

    // for (auto& v : _initial_vertices)
    //     v = rot_mat * v;

    _mesh_origin = rot_mat * _mesh_origin;
    moveTogether(p);
}

std::tuple<Real, Vec3r, Mat3r> ParticleMesh::massProperties(Real density) const
{
    // uses the algorithm described here: http://number-none.com/blow/inertia/index.html
    Real total_volume = 0;
    Vec3r weighted_volume(0,0,0);
    Mat3r covariance = Mat3r::Zero();

    // covariance of "canonical" tetrahedron which is (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    Mat3r C_canonical;
    C_canonical <<  1.0/60.0, 1.0/120.0, 1.0/120.0,
                    1.0/120.0, 1.0/60.0, 1.0/120.0,
                    1.0/120.0, 1.0/120.0, 1.0/60.0;
    for (const auto& f : _faces)
    {
        // each triangle in the mesh + origin forms a tetrahedron
        // v0=origin, v1=f[0], v2=f[1], v3=f[2]
        const Vec3r v0(0,0,0);
        const Vec3r v1 = vertex(f[0]);
        const Vec3r v2 = vertex(f[1]);
        const Vec3r v3 = vertex(f[2]);

        // tet basis matrix
        Mat3r A;
        A.col(0) = (v1 - v0);
        A.col(1) = (v2 - v0);
        A.col(2) = (v3 - v0);

        // find signed volume of tet
        const Real volume = A.determinant() / 6.0;

        // calculate the center of mass of this tetrahedron - just average of 4 vertices
        const Vec3r tet_cm = 0.25*(v0 + v1 + v2 + v3);
        // update overall center of mass using a weighted average
        weighted_volume += tet_cm*volume;

        // add covariance matrix from this tet
        covariance += A.determinant() * A * C_canonical * A.transpose();
        // update overall volume
        total_volume += volume;
    }

    Vec3r center_of_mass = weighted_volume / total_volume;

    // move covariance matrix to center of mass
    covariance = covariance + total_volume * ( 2*(-center_of_mass) * (center_of_mass).transpose() + (-center_of_mass)*(-center_of_mass).transpose());

    // compute moment of inertia tensor from covariance mat
    const Mat3r I = Mat3r::Identity() * covariance.trace() - covariance;

    return std::tuple<Real, Vec3r, Mat3r>(density*total_volume, center_of_mass, density*I);
}

Vec3r ParticleMesh::massCenter() const
{
    Real total_volume = 0;
    Vec3r weighted_volume(0,0,0);
    for (const auto& f : _faces)
    {
        // each triangle in the mesh + origin forms a tetrahedron
        // v0=origin, v1=f[0], v2=f[1], v3=f[2]
        const Vec3r v0(0,0,0);
        const Vec3r v1 = vertex(f[0]);
        const Vec3r v2 = vertex(f[1]);
        const Vec3r v3 = vertex(f[2]);

        // tet basis matrix
        Mat3r A;
        A.col(0) = (v1 - v0);
        A.col(1) = (v2 - v0);
        A.col(2) = (v3 - v0);

        // find signed volume of tet
        const Real volume = A.determinant() / 6.0;

        // calculate the center of mass of this tetrahedron - just average of 4 vertices
        const Vec3r tet_cm = 0.25*(v0 + v1 + v2 + v3);
        // update overall center of mass using a weighted average
        // if (total_volume + volume > 0)
        //     center_of_mass = (center_of_mass*total_volume + tet_cm*volume) / (total_volume + volume);
        weighted_volume += tet_cm * volume;

        // update overall volume
        total_volume += volume;
    }

    return weighted_volume / total_volume;
}

bool ParticleMesh::isInside(const Vec3r& p) const
{
    Real total = 0;
    for (const auto& f : _faces)
    {
        // compute solid angle for each face
        Vec3r a = vertex(f[0]) - p;
        Vec3r b = vertex(f[1]) - p;
        Vec3r c = vertex(f[2]) - p;

        Real la = a.norm();
        Real lb = b.norm();
        Real lc = c.norm();

        Real numerator = a.dot(b.cross(c));
        Real denominator = la*lb*lc + a.dot(b)*lc + b.dot(c)*la + c.dot(a)*lb;
        Real solid_angle = 2.0 * std::atan2(numerator, denominator);

        total += solid_angle;
    }
    // winding number = normalized total solid angle
    Real w = total / (4.0 * M_PI);

    // if winding number > 0.5, then the point is inside the mesh
    return std::abs(w) > 0.5;
}

void ParticleMesh::writeMeshToObjFile(const std::string& filename) const
{
    std::ofstream obj_file(filename);
    if (obj_file.is_open())
    {
        for (const auto& v : _vertices)
        {
            const Vec3r& p = _particle_pool->positions[v];
            obj_file << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
        }
        
        for (const auto& f : _faces)
        {
            obj_file << "f " << f[0]+1 << " " << f[1]+1 << " " << f[2]+1 << std::endl;
        }
    }
}