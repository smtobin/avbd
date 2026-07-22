#pragma once

#include "common/common.hpp"

/** Enumerates the different possible collision types in the sim */
enum class DetectedCollisionType : uint8_t
{
    TriangleTriangle_VertexFace,
    TriangleTriangle_EdgeEdge,
    TriangleRigid,
    RigidRigid
};

/** Type-tagged union representing a detected collision during narrow-phase collision detection.
 * Contains a unique key which is generated from the type of collision, the primitive IDs in collision, and the features involved in the collision.
 */
struct DetectedCollision
{
    DetectedCollisionType type; // type of collision
    uint64_t key;           // unique key identifying the detected collision
    unsigned gen1, gen2;    // generations of each of the primitives in the collision pool
    Vec3r normal;   // collision normal
    union {
        struct { unsigned tri_vertex; Vec3u tri; Vec3r barys;  } TriangleTriangle_VertexFace;
        struct { Vec2u edge1; Real s1; Vec2u edge2; Real s2; } TriangleTriangle_EdgeEdge;
        struct { Vec3u tri; Vec3r barys; unsigned rb; Vec3r cp_rb_local; } TriangleRigid;
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
        constexpr unsigned IDX_BITS  = 27;  // bits for primitive index: supports ~100 M primitives
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

};