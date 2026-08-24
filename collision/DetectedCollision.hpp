#pragma once

#include "common/common.hpp"

/** Enumerates the different possible collision types in the sim */
enum class DetectedCollisionType : uint8_t
{
    TriangleTriangle_VertexFace,
    TriangleTriangle_EdgeEdge,
    TriangleRigid,
    TriangleRod,
    RodRod,
    RodRigid,
    RigidRigid
};

/** Type-tagged union representing a detected collision during narrow-phase collision detection.
 * Contains a unique key which is generated from the type of collision, the primitive IDs in collision, and the features involved in the collision.
 */
struct DetectedCollision
{
    DetectedCollisionType type; // type of collision
    unsigned e_idx;         // the index of the collision energy in the appropriate energy pool (depending on the detected collision type)
    uint64_t key;           // unique key identifying the detected collision
    unsigned gen1, gen2;    // generations of each of the primitives in the collision pool
    Vec3r normal;   // collision normal
    union {
        struct { unsigned tri_vertex; Vec3u tri; Vec3r barys;  } TriangleTriangle_VertexFace;
        struct { Vec2u edge1; Real s1; Vec2u edge2; Real s2; } TriangleTriangle_EdgeEdge;
        struct { Vec3u tri; Vec3r barys; unsigned rb; Vec3r cp_rb_local; } TriangleRigid;
        struct { Vec3u tri; Vec3r barys; Vec2u rod; Real s; } TriangleRod;
        struct { Vec2u rod1; Real s1; Vec2u rod2; Real s2; } RodRod;
        struct { Vec2u rod; Real s; unsigned rb; Vec3r cp_rb_local; } RodRigid;
        struct { unsigned rb1; Vec3r cp_rb_local1; unsigned rb2; Vec3r cp_rb_local2; } RigidRigid;
    };

    /** Generate the key for a detected collision
     * @param type the type of detected collision
     * @param prim_id1 the index of the 1st primitive in collision
     * @param prim_id2 the index of the 2nd primitive in collision
     * @param code optional, the code corresponding to the specific features in collision
     */
    static inline uint64_t generateKey(DetectedCollisionType type, unsigned prim_idx1, unsigned prim_idx2, unsigned code=0)
    {
        // canonicalize the indices so that (A,B) and (B,A) hash the same
        unsigned idxA = prim_idx1, idxB = prim_idx2;
        if (idxA > idxB)
            std::swap(idxA, idxB);

        constexpr unsigned TYPE_BITS = 5;   // bits for collision type: supports up to 32 different collision types
        constexpr unsigned IDX_BITS  = 26;  // bits for primitive index: supports ~50 M primitives
        constexpr unsigned CODE_BITS = 7;   // bits for the collision code: supports 128 different codes
        static_assert(TYPE_BITS + 2 * IDX_BITS + CODE_BITS == 64, "key must fill 64 bits");

        constexpr unsigned IDX_MASK  = (1u << IDX_BITS)  - 1;
        constexpr unsigned CODE_MASK = (1u << CODE_BITS) - 1;

        assert(prim_idx1 <= IDX_MASK && prim_idx2 <= IDX_MASK && "primitive index overflow");
        assert(code <= CODE_MASK && "feature code overflow");

        uint64_t key = static_cast<uint64_t>(type) & ((1u << TYPE_BITS) - 1);
        key = (key << IDX_BITS)  | (idxA & IDX_MASK);
        key = (key << IDX_BITS)  | (idxB & IDX_MASK);
        key = (key << CODE_BITS) | (code & CODE_MASK);
        return key;
    }

    /** Default constructor */
    DetectedCollision() = default;

    /** Implement copy constructor and assignment */
private:
    void copyPayload(const DetectedCollision& other)
    {
        switch (other.type)
        {
        case DetectedCollisionType::TriangleTriangle_VertexFace:
            new (&TriangleTriangle_VertexFace)
                decltype(TriangleTriangle_VertexFace)(other.TriangleTriangle_VertexFace);
            break;

        case DetectedCollisionType::TriangleTriangle_EdgeEdge:
            new (&TriangleTriangle_EdgeEdge)
                decltype(TriangleTriangle_EdgeEdge)(other.TriangleTriangle_EdgeEdge);
            break;

        case DetectedCollisionType::TriangleRigid:
            new (&TriangleRigid)
                decltype(TriangleRigid)(other.TriangleRigid);
            break;

        case DetectedCollisionType::TriangleRod:
            new (&TriangleRod)
                decltype(TriangleRod)(other.TriangleRod);
            break;

        case DetectedCollisionType::RodRod:
            new (&RodRod)
                decltype(RodRod)(other.RodRod);
            break;
        
        case DetectedCollisionType::RodRigid:
            new (&RodRigid)
                decltype(RodRigid)(other.RodRigid);
            break;

        case DetectedCollisionType::RigidRigid:
            new (&RigidRigid)
                decltype(RigidRigid)(other.RigidRigid);
            break;
        }
    }
public:
    DetectedCollision(const DetectedCollision& other)
        : type(other.type),
        e_idx(other.e_idx),
        key(other.key),
        gen1(other.gen1),
        gen2(other.gen2),
        normal(other.normal)
    {
        copyPayload(other);
    }

    DetectedCollision& operator=(const DetectedCollision& other)
    {
        if (this == &other)
            return *this;

        type   = other.type;
        e_idx  = other.e_idx;
        key    = other.key;
        gen1   = other.gen1;
        gen2   = other.gen2;
        normal = other.normal;

        copyPayload(other);

        return *this;
    }

};