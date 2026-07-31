#pragma once

#include "common/common.hpp"
#include "common/WorkerThreadContext.hpp"
#include "common/ParticlePool.hpp"
#include "collision/CollisionPrimitivePool.hpp"
#include "collision/LBVH.hpp"

namespace Collision
{

struct LBVHBuilder
{
    static constexpr unsigned PARALLEL_RADIX_SORT_THRESHOLD = 50000;
    static constexpr unsigned PARALLEL_AABB_MORTON_THRESHOLD = 5000;
    static constexpr unsigned PARALLEL_CONSTRUCT_TREE_THRESHOLD = 500;
    static constexpr unsigned PARALLEL_ASSEMBLE_BVH_THRESHOLD = 5000;

    /** Scratch memory for parallel radix sort.
     * These will be resized appropriately within the function, they just need to exist outside the function.
     */
    std::vector<unsigned> _radix_sort_combined_offsets;
    std::vector<unsigned> _radix_sort_temp;

    void buildBVH(const ParticlePool& particle_pool, CollisionPrimitivePool& col_pool, LBVH& lbvh, Real dt);

    void computeAABB_MortonCode(const ParticlePool& particle_pool, CollisionPrimitivePool& col_pool, Real dt);
    void constructTree(CollisionPrimitivePool& col_pool, LBVH& lbvh);
    void assembleBVH(CollisionPrimitivePool& col_pool, LBVH& lbvh);

    void buildBVH_Parallel(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, Sim::SimulationContext& ctx);
    void computeAABB_MortonCode_Parallel(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, Sim::SimulationContext& ctx);
    void constructTree_Parallel(WorkerThreadContext& w_ctx, Sim::SimulationContext& ctx);
    void assembleBVH_Parallel(WorkerThreadContext& w_ctx, Sim::SimulationContext& ctx);

private:
    static inline int commonPrefixLen(uint64_t code_i, unsigned i, uint64_t code_j, unsigned j)
    {
        if (code_i == code_j)
        {
            // ties broken by index — extend LCP count by 64 + clz of index XOR
            return 64 + std::countl_zero(static_cast<uint64_t>(i) ^ static_cast<uint64_t>(j));
        }
        return std::countl_zero(code_i ^ code_j);
    }

    /** Morton code generation */
    static inline uint32_t expandBits10(uint32_t v)
    {
        v &= 0x3FFu;
        v = (v | (v << 16)) & 0x030000FFu;
        v = (v | (v << 8))  & 0x0300F00Fu;
        v = (v | (v << 4))  & 0x030C30C3u;
        v = (v | (v << 2))  & 0x09249249u;
        return v;
    }

    static inline uint32_t morton3D_32(const Vec3r& p_normalized)
    {
        Real x = std::clamp(p_normalized[0] * 1024.0, 0.0, 1023.0);
        Real y = std::clamp(p_normalized[1] * 1024.0, 0.0, 1023.0);
        Real z = std::clamp(p_normalized[2]* 1024.0, 0.0, 1023.0);
        uint32_t xx = expandBits10(static_cast<uint32_t>(x));
        uint32_t yy = expandBits10(static_cast<uint32_t>(y));
        uint32_t zz = expandBits10(static_cast<uint32_t>(z));
        return xx | (yy << 1) | (zz << 2);
    }

    static inline uint64_t expandBits21(uint64_t v)
    {
        v &= 0x1FFFFFull;
        v = (v | (v << 32)) & 0x1F00000000FFFFull;
        v = (v | (v << 16)) & 0x1F0000FF0000FFull;
        v = (v | (v << 8))  & 0x100F00F00F00F00Full;
        v = (v | (v << 4))  & 0x10C30C30C30C30C3ull;
        v = (v | (v << 2))  & 0x1249249249249249ull;
        return v;
    }

    static inline uint64_t morton3D_64(const Vec3r& p_normalized)
    {
        Real x = std::clamp(p_normalized[0] * 2097152.0, 0.0, 2097151.0); // 2^21
        Real y = std::clamp(p_normalized[1] * 2097152.0, 0.0, 2097151.0);
        Real z = std::clamp(p_normalized[2] * 2097152.0, 0.0, 2097151.0);
        uint64_t xx = expandBits21(static_cast<uint64_t>(x));
        uint64_t yy = expandBits21(static_cast<uint64_t>(y));
        uint64_t zz = expandBits21(static_cast<uint64_t>(z));
        return xx | (yy << 1) | (zz << 2);
    }

};

} // namespace Collision