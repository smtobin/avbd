#include "common/common.hpp"
#include "collision/LBVH.hpp"
#include "collision/LBVHBuilder.hpp"
#include "collision/LBVHTraversal.hpp"

int main()
{
    int capacity = 1000;
    Collision::CollisionPrimitivePool col_pool(capacity);

    Collision::LBVH bvh;
    /** Test radix tree construction using the example from Karras 2012 */
    std::array<uint64_t, 8> morton_codes = {
        0b11001ULL,
        0b00100ULL,
        0b00101ULL,
        0b00010ULL,
        0b11000ULL,
        0b00001ULL,
        0b10011ULL,
        0b11110ULL
    };

    for (unsigned i = 0; i < morton_codes.size(); i++)
    {
        unsigned new_slot = col_pool.allocSlot();
        std::cout << "New slot: " << new_slot << std::endl;
        col_pool.morton_code[new_slot] = morton_codes[i];
    }

    Collision::LBVHBuilder::radixSort(col_pool.morton_code, col_pool.sorted_order, col_pool.totalSize());
    std::cout << "Sorted order:" << std::endl;
    for (unsigned i = 0; i < col_pool.totalSize(); i++)
    {
        std::cout << " " << i << ": " << col_pool.sorted_order[i] << std::endl;
    }
    Collision::LBVHBuilder::constructTree(col_pool, bvh);
    std::cout << "Root: " << bvh.root << std::endl;
    for (unsigned i = 0; i < bvh.parent.size(); i++)
    {
        std::cout << "Node " << i << ": Left=" << bvh.left[i] << " Right=" << bvh.right[i] << " Parent=" << bvh.parent[i] << " Leaf count=" << bvh.leaf_count[i] << std::endl;
    }
    bvh.printTreeWithInfo(col_pool);

    


}