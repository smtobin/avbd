#include "collision/LBVHBuilder.hpp"
#include "common/Algorithm.hpp"
#include "simulation/SimulationContext.hpp"

namespace Collision
{

static constexpr unsigned PARALLEL_RADIX_SORT_THRESHOLD = 50000;
static constexpr unsigned PARALLEL_AABB_MORTON_THRESHOLD = 500;
static constexpr unsigned PARALLEL_CONSTRUCT_TREE_THRESHOLD = 500;
static constexpr unsigned PARALLEL_ASSEMBLE_BVH_THRESHOLD = 5000;

void LBVHBuilder::buildBVH(const ParticlePool& particle_pool, CollisionPrimitivePool& col_pool, LBVH& lbvh, Real dt)
{
    // std::cout << "Building BVH..." << std::endl;
    // compute AABBs and Morton codes
    computeAABB_MortonCode(particle_pool, col_pool, dt);
    // radix sort by Morton code
    Algorithm::radixSort(col_pool.morton_code, col_pool.sorted_order, col_pool.totalSize());
    // construct the radix tree
    constructTree(col_pool, lbvh);
    // coalesce bounding boxes using radix tree
    assembleBVH(col_pool, lbvh);
    // std::cout << "Done." << std::endl;
}

void LBVHBuilder::buildBVH_Parallel(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, Sim::SimulationContext& ctx)
{
    if (ctx.collision_pool.totalSize() < PARALLEL_AABB_MORTON_THRESHOLD)
    {
        if (w_ctx.idx == 0)
            computeAABB_MortonCode(ctx.particles, ctx.collision_pool, ctx.params.dt);
        w_ctx.barrier->arrive_and_wait();
    }
    else
    {
        computeAABB_MortonCode_Parallel(w_ctx, all_worker_contexts, ctx);

        if (w_ctx.idx == 0)
        {
            // auto morton_copy = ctx.collision_pool.morton_code;
            computeAABB_MortonCode(ctx.particles, ctx.collision_pool, ctx.params.dt);
            // for (unsigned i = 0; i < ctx.collision_pool.totalSize(); i++)
            // {
            //     std::cout << "Parallel Morton " << i << ": " << morton_copy[i] << " Serial: " << ctx.collision_pool.morton_code[i] << std::endl;
            // }
        }
            
        w_ctx.barrier->arrive_and_wait();
    }

    if (ctx.collision_pool.totalSize() < PARALLEL_RADIX_SORT_THRESHOLD)
    {
        if (w_ctx.idx == 0)
            Algorithm::radixSort(ctx.collision_pool.morton_code, ctx.collision_pool.sorted_order, ctx.collision_pool.totalSize());
        w_ctx.barrier->arrive_and_wait();
    }
    else
    {
        Algorithm::radixSort_Parallel(
            w_ctx, 
            all_worker_contexts, 
            _radix_sort_combined_offsets, 
            ctx.collision_pool.morton_code, 
            ctx.collision_pool.sorted_order, 
            _radix_sort_temp, 
            ctx.collision_pool.totalSize()
        );
    }

    if (ctx.collision_pool.totalSize() < PARALLEL_CONSTRUCT_TREE_THRESHOLD)
    {
        if (w_ctx.idx == 0)
            constructTree(ctx.collision_pool, ctx.lbvh);
        w_ctx.barrier->arrive_and_wait();
    }
    else
    {
        constructTree_Parallel(w_ctx, ctx);
    }

    if (ctx.collision_pool.totalSize() < PARALLEL_ASSEMBLE_BVH_THRESHOLD)
    {
        if (w_ctx.idx == 0)
            assembleBVH(ctx.collision_pool, ctx.lbvh);
        w_ctx.barrier->arrive_and_wait();
    }
    else
    {
        assembleBVH_Parallel(w_ctx, ctx);
    }

    

    
}

void LBVHBuilder::computeAABB_MortonCode(const ParticlePool& particle_pool, CollisionPrimitivePool& col_pool, Real dt)
{
    // iterate through primtivies and compute AABB, centroid, and Morton code
    AABB scene_box = AABB::empty();
    for (unsigned p_idx : col_pool)
    {
        col_pool.aabb[p_idx] = col_pool.speculativeGlobalBounds(p_idx, particle_pool, dt);

        // centroid
        col_pool.centroid[p_idx] = col_pool.aabb[p_idx].center();

        scene_box.expand(col_pool.aabb[p_idx]);

        // std::cout << "Centroid for primitive " << p_idx << ": " << col_pool.centroid[p_idx].transpose() << std::endl;
        // std::cout << "AABB for primitive " << p_idx << ": " << col_pool.aabb[p_idx] << std::endl;
    }

    // Morton code
    Vec3r extent = scene_box.max - scene_box.min;
    extent = extent.cwiseMax(1e-6); // guard against degenerate axes
    for (unsigned p_idx : col_pool)
    {
        // normalize centroid between [0,1]^3
        Vec3r normalized_centroid = (col_pool.centroid[p_idx] - scene_box.min).cwiseQuotient(extent);
        col_pool.morton_code[p_idx] = morton3D_64(normalized_centroid);

        // std::cout << "Normalized centroid for primitive " << p_idx << ": " << normalized_centroid.transpose() << std::endl;
        // std::cout << "Morton code for primitive " << p_idx << ": " << col_pool.morton_code[p_idx] << std::endl;
    }
}

void LBVHBuilder::computeAABB_MortonCode_Parallel(WorkerThreadContext& w_ctx, const std::vector<WorkerThreadContext>& all_worker_contexts, Sim::SimulationContext& ctx)
{
    // std::cout << "Num primitives: " << ctx.collision_pool.totalSize() << std::endl;
    auto [start, end] = w_ctx.computeStartEnd(ctx.collision_pool);

    // reset this thread's scene box
    w_ctx.BuildBVHContext.scene_box = AABB::empty();

    for (unsigned p_idx = start; p_idx < end; p_idx++)
    {
        if (!ctx.collision_pool.active[p_idx])
            continue;

        ctx.collision_pool.aabb[p_idx] = ctx.collision_pool.speculativeGlobalBounds(p_idx, ctx.particles, ctx.params.dt);

        // centroid
        ctx.collision_pool.centroid[p_idx] = ctx.collision_pool.aabb[p_idx].center();   /** TODO: (07/30/26) change to centroid of primitive? */

        w_ctx.BuildBVHContext.scene_box.expand(ctx.collision_pool.aabb[p_idx]);
    }

    // wait for all threads to complete and coalesce thread-local scene boxes
    w_ctx.barrier->arrive_and_wait();

    // combine local scene boxes
    if (w_ctx.idx == 0)
    {
        ctx.scene_box = w_ctx.BuildBVHContext.scene_box;
        for (unsigned other_idx = 1; other_idx < all_worker_contexts.size(); other_idx++)
        {
            ctx.scene_box.expand(all_worker_contexts[other_idx].BuildBVHContext.scene_box);
        }
    }

    w_ctx.barrier->arrive_and_wait();

    // compute Morton codes
    Vec3r extent = ctx.scene_box.max - ctx.scene_box.min;
    extent = extent.cwiseMax(1e-6); // guard against degenerate axes
    for (unsigned p_idx = start; p_idx < end; p_idx++)
    {
        if (!ctx.collision_pool.active[p_idx])
            continue;

        // normalize centroid between [0,1]^3
        Vec3r normalized_centroid = (ctx.collision_pool.centroid[p_idx] - ctx.scene_box.min).cwiseQuotient(extent);
        ctx.collision_pool.morton_code[p_idx] = morton3D_64(normalized_centroid);
    }

    w_ctx.barrier->arrive_and_wait();
}

void LBVHBuilder::constructTree(CollisionPrimitivePool& col_pool, LBVH& lbvh)
{
    int n = static_cast<int>(col_pool.totalSize());
    lbvh.resize(n);
    lbvh.parent.assign(lbvh.parent.size(), LBVH::INVALID);
    lbvh.left.assign(lbvh.left.size(), LBVH::INVALID);
    lbvh.right.assign(lbvh.right.size(), LBVH::INVALID);

    auto delta = [&](int i, int j)
    {
        if (j < 0 || j >= n)
            return -1;

        return commonPrefixLen(
            col_pool.morton_code[col_pool.sorted_order[i]],
            col_pool.sorted_order[i], 
            col_pool.morton_code[col_pool.sorted_order[j]],
            col_pool.sorted_order[j]
        );
    };

    // iterate over internal nodes
    for (int i = 0; i < n-1; i++)
    {
        // "direction" of interval
        int d = (delta(i, i + 1) - delta(i, i - 1) >= 0) ? 1 : -1;
        
        // lower bound on LCP length for siblings 
        int dmin = delta(i, i-d);

        // find upper bound on the range length
        int l_max = 128;
        while (delta(i, i+l_max*d) > dmin)
        {
            l_max *= 4;
        }

        // binary search for the exact far end
        int l = 0;
        for (unsigned t = l_max; t >= 1; t/=2)
        {
            if (delta(i, i + (l+t)*d) > dmin)
                l += t;
        }

        int j = i + l*d;   // range end

        // binary search for the split position within [i, j]
        int dnode = delta(i, j);
        int s = 0;
        int div = 2;
        for (int t = (l + div - 1) / div; t >= 1; t = (t == 1 ? 0 : (t + div - 1) / div))
        {
            int new_s = s + t;
            if (new_s < l && delta(i, i + new_s * d) > dnode)
                s = new_s;
        }
        int gamma = i + s * d + std::min(d, 0);


        // index of the current internal node in the global LBVH arrays
        unsigned global_idx = static_cast<unsigned>(i) + n;

        if (std::min(i, j) == gamma)
        {
            lbvh.left[global_idx] = static_cast<unsigned>(gamma);
            lbvh.parent[gamma] = global_idx;
        }
        else
        {
            lbvh.left[global_idx] = static_cast<unsigned>(n + gamma);
            lbvh.parent[n+gamma] = global_idx;
        }

        if (std::max(i, j) == gamma + 1)
        {
            lbvh.right[global_idx] = static_cast<unsigned>(gamma + 1);
            lbvh.parent[gamma + 1] = global_idx;
        }
        else
        {
            lbvh.right[global_idx] = static_cast<unsigned>(n + gamma + 1);
            lbvh.parent[n + gamma + 1] = global_idx;
        }
    }

    // find root - start at an arbitrary leaf node and go upwards
    unsigned node = 0;
    while (lbvh.parent[node] != LBVH::INVALID)
        node = lbvh.parent[node];
    lbvh.root = node;

    // process leaves
    for (int l_idx = 0; l_idx < n; l_idx++)
    {
        lbvh.leaf_start[l_idx] = l_idx;
        lbvh.leaf_count[l_idx] = 1;
        lbvh.subtree_size[l_idx] = 1;
    }

    // internal nodes
    for (int idx = 0; idx+1 < n; idx++)
    {
        lbvh.leaf_count[n + idx] = 0;
    }
}

void LBVHBuilder::constructTree_Parallel(WorkerThreadContext& w_ctx, Sim::SimulationContext& ctx)
{
    int n = static_cast<int>(ctx.collision_pool.totalSize());
    LBVH& lbvh = ctx.lbvh;

    // allocate memory and reset the BVH
    if (w_ctx.idx == 0)
    {
        lbvh.resize(n);
        lbvh.parent.assign(lbvh.parent.size(), LBVH::INVALID);
        lbvh.left.assign(lbvh.left.size(), LBVH::INVALID);
        lbvh.right.assign(lbvh.right.size(), LBVH::INVALID);
    }
    w_ctx.barrier->arrive_and_wait();

    auto delta = [&](int i, int j)
    {
        if (j < 0 || j >= n)
            return -1;

        return commonPrefixLen(
            ctx.collision_pool.morton_code[ctx.collision_pool.sorted_order[i]],
            ctx.collision_pool.sorted_order[i], 
            ctx.collision_pool.morton_code[ctx.collision_pool.sorted_order[j]],
            ctx.collision_pool.sorted_order[j]
        );
    };

    // iterate over internal nodes
    auto [istart, iend] = w_ctx.computeStartEnd(n-1);
    for (int i = istart; i < static_cast<int>(iend); i++)
    {
        // "direction" of interval
        int d = (delta(i, i + 1) - delta(i, i - 1) >= 0) ? 1 : -1;
        
        // lower bound on LCP length for siblings 
        int dmin = delta(i, i-d);

        // find upper bound on the range length
        int l_max = 128;
        while (delta(i, i+l_max*d) > dmin)
        {
            l_max *= 4;
        }

        // binary search for the exact far end
        int l = 0;
        for (unsigned t = l_max; t >= 1; t/=2)
        {
            if (delta(i, i + (l+t)*d) > dmin)
                l += t;
        }

        int j = i + l*d;   // range end

        // binary search for the split position within [i, j]
        int dnode = delta(i, j);
        int s = 0;
        int div = 2;
        for (int t = (l + div - 1) / div; t >= 1; t = (t == 1 ? 0 : (t + div - 1) / div))
        {
            int new_s = s + t;
            if (new_s < l && delta(i, i + new_s * d) > dnode)
                s = new_s;
        }
        int gamma = i + s * d + std::min(d, 0);


        // index of the current internal node in the global LBVH arrays
        unsigned global_idx = static_cast<unsigned>(i) + n;

        if (std::min(i, j) == gamma)
        {
            lbvh.left[global_idx] = static_cast<unsigned>(gamma);
            lbvh.parent[gamma] = global_idx;
        }
        else
        {
            lbvh.left[global_idx] = static_cast<unsigned>(n + gamma);
            lbvh.parent[n+gamma] = global_idx;
        }

        if (std::max(i, j) == gamma + 1)
        {
            lbvh.right[global_idx] = static_cast<unsigned>(gamma + 1);
            lbvh.parent[gamma + 1] = global_idx;
        }
        else
        {
            lbvh.right[global_idx] = static_cast<unsigned>(n + gamma + 1);
            lbvh.parent[n + gamma + 1] = global_idx;
        }
    }

    // process leaves
    auto [lstart, lend] = w_ctx.computeStartEnd(n);
    for (int l_idx = lstart; l_idx < static_cast<int>(lend); l_idx++)
    {
        lbvh.leaf_start[l_idx] = l_idx;
        lbvh.leaf_count[l_idx] = 1;
        lbvh.subtree_size[l_idx] = 1;
    }

    // internal nodes
    for (int idx = istart; idx < static_cast<int>(iend); idx++)
    {
        lbvh.leaf_count[n + idx] = 0;
    }

    w_ctx.barrier->arrive_and_wait();
    
    // let thread 0 find root - start at an arbitrary leaf node and go upwards
    if (w_ctx.idx == 0)
    {
        unsigned node = 0;
        while (lbvh.parent[node] != LBVH::INVALID)
            node = lbvh.parent[node];
        lbvh.root = node;
    }

}

void LBVHBuilder::assembleBVH(CollisionPrimitivePool& col_pool, LBVH& lbvh)
{
    unsigned n = col_pool.totalSize();

    // refit pass
    std::vector<std::atomic<uint8_t>> visited(2 * n - 1);

    // start at leaves, then walk up
    for (unsigned l_idx = 0; l_idx < n; l_idx++)
    {
        AABB box = AABB::empty();
        for (unsigned i = 0; i < lbvh.leaf_count[l_idx]; i++)
        {
            unsigned prim = col_pool.sorted_order[lbvh.leaf_start[l_idx] + i];
            box.expand(col_pool.aabb[prim]);
        }

        lbvh.min_x[l_idx] = box.min[0]; lbvh.min_y[l_idx] = box.min[1]; lbvh.min_z[l_idx] = box.min[2];
        lbvh.max_x[l_idx] = box.max[0]; lbvh.max_y[l_idx] = box.max[1]; lbvh.max_z[l_idx] = box.max[2];

        unsigned node = lbvh.parent[l_idx];
        while (node != LBVH::INVALID)
        {
            if (visited[node].fetch_add(1) == 0)
                break; // first thread here - the sibling will finish the union

            // merge the AABB
            unsigned l = lbvh.left[node];
            unsigned r = lbvh.right[node];

            lbvh.min_x[node] = std::min(lbvh.min_x[l], lbvh.min_x[r]);
            lbvh.min_y[node] = std::min(lbvh.min_y[l], lbvh.min_y[r]);
            lbvh.min_z[node] = std::min(lbvh.min_z[l], lbvh.min_z[r]);
            lbvh.max_x[node] = std::max(lbvh.max_x[l], lbvh.max_x[r]);
            lbvh.max_y[node] = std::max(lbvh.max_y[l], lbvh.max_y[r]);
            lbvh.max_z[node] = std::max(lbvh.max_z[l], lbvh.max_z[r]);

            lbvh.subtree_size[node] = lbvh.subtree_size[l] + lbvh.subtree_size[r];

            node = lbvh.parent[node];
        }
        
    }
}

void LBVHBuilder::assembleBVH_Parallel(WorkerThreadContext& w_ctx, Sim::SimulationContext& ctx)
{
    /** TODO: (07/31/26) Do level-synchronous refit instead. Calculate depth of each node, then process each depth level in parallel. */
    unsigned n = ctx.collision_pool.totalSize();
    LBVH& lbvh = ctx.lbvh;

    // refit pass
    // visited flags were assigned 0 when BVH memory was allocated
    // start at the leaves, then walk up
    auto [lstart, lend] = w_ctx.computeStartEnd(n);
    for (unsigned l_idx = lstart; l_idx < lend; l_idx++)
    {
        // merge the bounding boxes of all the leaf primitives under this BVH leaf node
        AABB box = AABB::empty();
        for (unsigned i = 0; i < lbvh.leaf_count[l_idx]; i++)
        {
            unsigned prim = ctx.collision_pool.sorted_order[lbvh.leaf_start[l_idx] + i];
            box.expand(ctx.collision_pool.aabb[prim]);
        }
        lbvh.min_x[l_idx] = box.min[0]; lbvh.min_y[l_idx] = box.min[1]; lbvh.min_z[l_idx] = box.min[2];
        lbvh.max_x[l_idx] = box.max[0]; lbvh.max_y[l_idx] = box.max[1]; lbvh.max_z[l_idx] = box.max[2];

        // traverse up the tree, merging AABBs
        unsigned node = lbvh.parent[l_idx];
        while (node != LBVH::INVALID)
        {
            if (lbvh.visited[node].fetch_add(1) == 0)
                break;  // first thread here, the sibling will finish the union

            // merge the AABB
            unsigned l = lbvh.left[node];
            unsigned r = lbvh.right[node];

            lbvh.min_x[node] = std::min(lbvh.min_x[l], lbvh.min_x[r]);
            lbvh.min_y[node] = std::min(lbvh.min_y[l], lbvh.min_y[r]);
            lbvh.min_z[node] = std::min(lbvh.min_z[l], lbvh.min_z[r]);
            lbvh.max_x[node] = std::max(lbvh.max_x[l], lbvh.max_x[r]);
            lbvh.max_y[node] = std::max(lbvh.max_y[l], lbvh.max_y[r]);
            lbvh.max_z[node] = std::max(lbvh.max_z[l], lbvh.max_z[r]);

            lbvh.subtree_size[node] = lbvh.subtree_size[l] + lbvh.subtree_size[r];

            node = lbvh.parent[node];
        }
    }
}

} // namespace Collision