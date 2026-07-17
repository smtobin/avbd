#include "collision/LBVHBuilder.hpp"

namespace Collision
{

void LBVHBuilder::buildBVH(const ParticlePool& particle_pool, CollisionPrimitivePool& col_pool, LBVH& lbvh)
{
    // compute AABBs and Morton codes
    computeAABB_MortonCode(particle_pool, col_pool);

    // radix sort by Morton code
    radixSort(col_pool.morton_code, col_pool.sorted_order);

    // construct the radix tree
    constructTree(col_pool, lbvh);

    // coalesce bounding boxes using radix tree
    assembleBVH(col_pool, lbvh);
}

void LBVHBuilder::computeAABB_MortonCode(const ParticlePool& particle_pool, CollisionPrimitivePool& col_pool)
{
    // iterate through primtivies and compute AABB, centroid, and Morton code
    AABB scene_box;
    for (unsigned p_idx : col_pool)
    {
        // AABB
        AABB box;
        for (unsigned k = 0; k < col_pool.num_particles[p_idx]; k++)
        {
            box.expand(particle_pool.positions[col_pool.particle_indices[p_idx][k]]);
        }
        col_pool.aabb[p_idx] = box;

        // centroid
        col_pool.centroid[p_idx] = box.center();

        scene_box.expand(box);
    }

    // Morton code
    Vec3r extent = scene_bounds.max - scene_bounds.min;
    extent = extent.cwiseMax(Vec3r(1e-6)); // guard against degenerate axes
    for (unsigned p_idx : col_pool)
    {
        // normalize centroid between [0,1]^3
        Vec3r normalized_centroid = (col_pool.centroid[p_idx] - scene_box.min).cwiseQuotient(extent);
        // col_pool.morton_code[p_idx] = morton3D_32(normalized_centroid);
        col_pool.morton_code[p_idx] = morton3D_64(normalized_centroid);
    }
}

void LBVHBuilder::radixSort(const std::vector<uint64_t>& unsorted, std::vector<unsigned>& sorted_order)
{

    sorted_order.resize(unsorted.size());
    std::iota(sorted_order.begin(), sorted_order.end(), 0);

    std::vector<unsigned> temp(unsorted.size());

    constexpr int BITS = 8;
    constexpr int BUCKETS = 1 << BITS;
    constexpr int MASK = BUCKETS - 1;

    for (int shift = 0; shift < 64; shift += BITS)
    {
        size_t count[BUCKETS] = {};

        // count
        for (unsigned idx : sorted_order)
        {
            uint64_t x = unsorted[idx];
            ++count[(x >> shift) & MASK];
        }
        
        // prefix sums
        size_t sum = 0;
        for (int i = 0; i < BUCKETS; i++)
        {
            size_t c = count[i];
            count[i] = sum;
            sum += c;
        }

        // scatter
        for (unsigned idx : sorted_order)
        {
            uint64_t x = unsorted[idx];
            temp[count[(x >> shift) & MASK]++] = idx;
        }

        sorted_order.swap(temp);
    }
}

void LBVHBuilder::constructTree(CollisionPrimitivePool& col_pool, LBVH& lbvh)
{
    unsigned n = col_pool.morton_code.size();
    lbvh.resize(n);
    lbvh.parent.assign(lbvh.parent.size(), INVALID);

    auto delta = [&](unsigned i, unsigned j)
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
    for (unsigned i = 0; i < internal.size(); i++)
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
        unsigned global_idx = i+n;

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
}

void LBVHBuilder::assembleBVH(CollisionPrimitivePool& col_pool, LBVH& lbvh)
{
    // process leaves
    unsigned n = col_pool.morton_code.size();
    for (unsigned l_idx = 0; l_idx < n; l_idx++)
    {
        lbvh.leaf_start[l_idx] = l_idx;
        lbvh.leaf_count[l_idx] = 1;
        lbvh.subtree_size[l_idx] = 1;
    }

    // internal nodes
    for (unsigned idx = 0; idx < n-1; idx++)
    {
        lbvh.leaf_count[n + idx] = 0;
    }

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
        while (node != INVALID)
        {
            if (visited[node].fetch_add(1) == 0)
                return; // first thread here - the sibling will finish the union

            // merge the AABB
            unsigned l = lbvh.left[node];
            unsigned r = lbvh.right[node];

            lbvh.min_x[node] = std::min(lbvh.min_x[l], lbvh.min_x[r]);
            lbvh.min_y[node] = std::min(lbvh.min_y[l], lbvh.min_y[r]);
            lbvh.min_z[node] = std::min(lbvh.min_z[l], lbvh.min_z[r]);
            lbvh.max_x[node] = std::min(lbvh.max_x[l], lbvh.max_x[r]);
            lbvh.max_y[node] = std::min(lbvh.max_y[l], lbvh.max_y[r]);
            lbvh.max_z[node] = std::min(lbvh.max_z[l], lbvh.max_z[r]);

            lbvh.subtree_size[node] = lbvh.subtree_size[l] + lbvh.subtree_size[r];

            node = lbvh.parent[node];
        }
        
    }
}

} // namespace Collision