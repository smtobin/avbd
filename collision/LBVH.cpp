#include "collision/LBVH.hpp"
#include "collision/CollisionPrimitivePool.hpp"

#include <bitset>
#include <queue>

namespace Collision
{

void LBVH::_printTreeWithInfoImpl(unsigned node, const CollisionPrimitivePool& pool, const std::string& prefix, bool is_left) const
{
    std::string node_str;
    if (node >= numPrimitives())
        node_str = "I" + std::to_string(node - numPrimitives());
    else
        node_str = "L" + std::to_string(node);

    std::stringstream morton_code_ss;
    if (node < numPrimitives())
        morton_code_ss << "(code=" << std::bitset<64>(pool.morton_code[pool.sorted_order[node]]) << ")";

    std::cout << prefix
            << (is_left ? "├── " : "└── ")
            << node_str << " " << morton_code_ss.str() << 
            '\n';

    if (leaf_count[node] > 0) return;

    _printTreeWithInfoImpl(left[node], pool,
            prefix + (is_left ? "│   " : "    "),
            true);

    _printTreeWithInfoImpl(right[node], pool,
            prefix + (is_left ? "│   " : "    "),
            false);
}

std::vector<unsigned> LBVH::nodeDepths() const
{
    std::cout << "Parent size: " << parent.size() << std::endl;
    std::vector<unsigned> depth(parent.size(), 0);

    std::queue<unsigned> q;
    q.push(root);

    while (!q.empty())
    {
        unsigned n = q.front();
        q.pop();

        if (left[n] != INVALID)
        {
            depth[left[n]] = depth[n] + 1;
            q.push(left[n]);
        }

        if (right[n] != INVALID)
        {
            depth[right[n]] = depth[n] + 1;
            q.push(right[n]);
        }
    }

    return depth;
}

} // namespace Collision