#include "collision/LBVHBuilder.hpp"

namespace Collision
{

void LBVHBuilder::buildBVH(const ParticlePool& particle_pool, const OrientedParticlePool& oriented_particle_pool, CollisionPrimitivePool& col_pool, LBVH& lbvh)
{
    // std::cout << "Building BVH..." << std::endl;
    // compute AABBs and Morton codes
    computeAABB_MortonCode(particle_pool, oriented_particle_pool, col_pool);
    // radix sort by Morton code
    radixSort(col_pool.morton_code, col_pool.sorted_order, col_pool.totalSize());
    // construct the radix tree
    constructTree(col_pool, lbvh);
    // coalesce bounding boxes using radix tree
    assembleBVH(col_pool, lbvh);
    // std::cout << "Done." << std::endl;
}

void LBVHBuilder::computeAABB_MortonCode(const ParticlePool& particle_pool, const OrientedParticlePool& oriented_particle_pool, CollisionPrimitivePool& col_pool)
{
    // iterate through primtivies and compute AABB, centroid, and Morton code
    AABB scene_box = AABB::empty();
    for (unsigned p_idx : col_pool)
    {
        col_pool.aabb[p_idx] = col_pool.globalBounds(p_idx, particle_pool, oriented_particle_pool);

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

void LBVHBuilder::radixSort(const std::vector<uint64_t>& unsorted, std::vector<unsigned>& sorted_order, unsigned size)
{

    sorted_order.resize(size);
    std::iota(sorted_order.begin(), sorted_order.end(), 0);

    std::vector<unsigned> temp(size);

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

} // namespace Collision