#pragma once

#include "common/common.hpp"
#include "common/Math.hpp"

#include "common/ParticlePool.hpp"

namespace Collision
{

struct TriangleTriangleCollider
{

    struct TriangleTriangleContact
    {
        bool edge_edge;
        union 
        {
            struct { unsigned a1,a2,b1,b2; Real sa,sb; } EdgeEdgeContact;
            struct { unsigned a1,b1,b2,b3; } PointFaceContact;
        };
    };

    struct Interval
    {
        Real min, max;
    };

    /** Narrow-phase collision detection between two triangles.
     * @param pool the particle memory pool that stores the vertices
     * @param a_idx1,a_idx2,a_idx3 the particle indices for triangle A
     * @param b_idx1,b_idx2,b_idx3 the particle indices for triangle B
     * @param tri_contact (output) the triangle-triangle contact
     * @returns true if there is contact, false otherwise
     */
    static bool collideTriangles(
        const ParticlePool& pool,
        unsigned a_idx1, unsigned a_idx2, unsigned a_idx3,
        unsigned b_idx1, unsigned b_idx2, unsigned b_idx3,
        TriangleTriangleContact& tri_contact
    )
    {
        /** TODO: (07/20/26) check for coplanarity and test this subroutine */

        // extract vertices from particle pool
        // triangle A
        const Vec3r& a1 = pool.positions[a_idx1];
        const Vec3r& a2 = pool.positions[a_idx2];
        const Vec3r& a3 = pool.positions[a_idx3];
        // triangle B
        const Vec3r& b1 = pool.positions[b_idx1];
        const Vec3r& b2 = pool.positions[b_idx2];
        const Vec3r& b3 = pool.positions[b_idx3];

        // compute centroids of each tri
        Vec3r c_a = (a1 + a2 + a3)/3.0;
        Vec3r c_b = (b1 + b2 + b3)/3.0;
        
        // separating axis theorem - look for axes where projections don't overlap
        // 11 potential separating axes:
        //  - 2 face normals
        //  - 9 edge x edge

        Real min_penetration = std::numeric_limits<Real>::max();
        Vec3r normal;
        bool flip_normal;
        int code;

        auto test_axis = [&](const Vec3r& axis, int cur_code) -> bool
        {
            // project triangle A onto axis
            Interval proj_a = _projectTriangle(a1, a2, a3, axis);
            // project triangle B onto axis
            Interval proj_b = _projectTriangle(b1, b2, b3, axis); 

            // find overlap between projections
            Real overlap = std::min(proj_a.max, proj_b.max) - std::max(proj_a.min, proj_b.min);

            if (overlap < 0)
                return false;   // separating axis found
            
            if (overlap < min_penetration)
            {
                min_penetration = overlap;
                normal = axis;
                code = cur_code;
                flip_normal = (c_b - c_a).dot(axis) < 0;
            }

            return true;
        };

        auto test_edge_edge_axis = [&](const Vec3r& e1, const Vec3r& e2, int cur_code) -> bool
        {
            Vec3r axis = e1.cross(e2);
            Real sq_len = axis.squaredNorm();
            // skip near-parallel edges
            if (sq_len < 1e-6)
                return true;

            axis = axis / std::sqrt(sq_len);

            return test_axis(axis, cur_code);
        };

        // triangle edges
        Vec3r e12_a = a2 - a1; Vec3r e13_a = a3 - a1; Vec3r e23_a = a3 - a2;
        Vec3r e12_b = b2 - b1; Vec3r e13_b = b3 - b1; Vec3r e23_b = b3 - b2;
        
        // test edge cross product axes with SAT
        if (!test_edge_edge_axis(e12_a, e12_b, 0)) return false;
        if (!test_edge_edge_axis(e12_a, e13_b, 1)) return false;
        if (!test_edge_edge_axis(e12_a, e23_b, 2)) return false;
        if (!test_edge_edge_axis(e13_a, e12_b, 3)) return false;
        if (!test_edge_edge_axis(e13_a, e13_b, 4)) return false;
        if (!test_edge_edge_axis(e13_a, e23_b, 5)) return false;
        if (!test_edge_edge_axis(e23_a, e12_b, 6)) return false;
        if (!test_edge_edge_axis(e23_a, e13_b, 7)) return false;
        if (!test_edge_edge_axis(e23_a, e23_b, 8)) return false;

        // test face axes with SAT
        Vec3r n_a = e12_a.cross(e13_a);
        n_a /= n_a.norm();
        Vec3r n_b = e12_b.cross(e13_b);
        n_b /= n_b.norm();

        if (!test_axis(n_a, 9)) return false;
        if (!test_axis(n_b, 10)) return false;

        // orient the normal so that it always points from triangle 1 -> triangle 2
        if (flip_normal)
            normal = -normal;

        // edge-edge contact
        if (code < 9)
        {
            Vec3r pa1, pa2, pb1, pb2;   // vertices that make up edges
            unsigned pai1, pai2, pbi1, pbi2;    // indices that make up edges
            int code_a = code / 3;
            int code_b = code % 3;
            switch(code_a)
            {
                case 0:
                    pa1 = a1;
                    pa2 = a2;
                    pai1 = a_idx1;
                    pai2 = a_idx2;
                    break;
                case 1:
                    pa1 = a1;
                    pa2 = a3;
                    pai1 = a_idx1;
                    pai2 = a_idx3;
                    break;
                case 2:
                    pa1 = a2;
                    pa2 = a3;
                    pai1 = a_idx2;
                    pai2 = a_idx3;
                    break;
            }
            switch (code_b)
            {
                case 0:
                    pb1 = b1;
                    pb2 = b2;
                    pbi1 = b_idx1;
                    pbi2 = b_idx2;
                    break;
                case 1:
                    pb1 = b1;
                    pb2 = b3;
                    pbi1 = b_idx1;
                    pbi2 = b_idx3;
                    break;
                case 2:
                    pb1 = b2;
                    pb2 = b3;
                    pbi1 = b_idx2;
                    pbi2 = b_idx3;
                    break;
            }

            // compute closest points between two edges
            Real sa, sb;
            Math::closestPoint_SegmentSegment(pa1, pa2, pb1, pb2, sa, sb);
            
            // create the contact
            tri_contact.edge_edge = true;
            tri_contact.EdgeEdgeContact = { pai1, pai2, pbi1, pbi2, sa, sb };
            return true;
        }
        // face contact
        else
        {
            // code 9 = triangle A's normal is MSA
            // code 10 = triangle B's normal is MSA
            bool a_is_ref = (code == 9);

            // reference triangle (the face)
            const Vec3r& r1  = a_is_ref ? a1 : b1;
            unsigned     ri1 = a_is_ref ? a_idx1 : b_idx1;
            unsigned     ri2 = a_is_ref ? a_idx2 : b_idx2;
            unsigned     ri3 = a_is_ref ? a_idx3 : b_idx3;

            // incident triangle (contributes the point)
            const Vec3r& i1  = a_is_ref ? b1 : a1;
            const Vec3r& i2  = a_is_ref ? b2 : a2;
            const Vec3r& i3  = a_is_ref ? b3 : a3;
            unsigned     ii1 = a_is_ref ? b_idx1 : a_idx1;
            unsigned     ii2 = a_is_ref ? b_idx2 : a_idx2;
            unsigned     ii3 = a_is_ref ? b_idx3 : a_idx3;

            // normal already points from the reference face toward the incident triangle
            // so the deepest-penetrating incident vertex is just the one with the smallest projection onto it
            Real d1 = (i1 - r1).dot(normal);
            Real d2 = (i2 - r1).dot(normal);
            Real d3 = (i3 - r1).dot(normal);

            unsigned deep_idx = ii1;
            Real deep_d = d1;
            if (d2 < deep_d) { deep_d = d2; deep_idx = ii2; }
            if (d3 < deep_d) { deep_d = d3; deep_idx = ii3; }

            tri_contact.edge_edge = false;
            tri_contact.PointFaceContact = { deep_idx, ri1, ri2, ri3 };
            return true;
        }

    }

private:

    /** Project triangle onto axis */
    static inline Interval _projectTriangle(
        const Vec3r& v1, const Vec3r& v2, const Vec3r& v3,
        const Vec3r& axis
    )
    {
        Real p1 = v1.dot(axis);
        Real p2 = v2.dot(axis);
        Real p3 = v3.dot(axis);

        return { std::min({p1, p2, p3}), std::max({p1, p2, p3}) };
    }

};

} // namespace Collision