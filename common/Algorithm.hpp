#pragma once

#include "common/common.hpp"

#include <numeric>

struct Algorithm
{

template<typename T>
static T lerp(const T& start, const T& end, Real interp)
{
    return (1-interp)*start + interp*end;
}

/** Radix sort on uint64 keys.
 * @param unsorted the list of unsorted keys
 * @param sorted_order (output) the indices of the unsorted array, arranged in the sorted order
 * @param size the size of the vector (the unsorted vector may have extra padding)
 */
static void radixSort(const std::vector<uint64_t>& unsorted, std::vector<unsigned>& sorted_order, unsigned size)
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
            ++count[(unsorted[idx] >> shift) & MASK];
        
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
            temp[count[(unsorted[idx] >> shift) & MASK]++] = idx;

        sorted_order.swap(temp);
    }
}

/** Radix sort on uint64 keys with a generic KeyFn for accessing the keys.
 * E.g. T is a struct type that has a "key" member variable, which KeyFn accesses.
 * @param data the list of unsorted structs
 * @param sorted_order (output) the indices of the unsorted array, arranged in the sorted order
 * @param size the size of the vector (unsorted vector may have extra padding)
 * @param key lambda function for extracting the key from T
 * 
 * Example usage:
 *      radixSort(collisions, sorted_order, collisions.size(), [](const DetectedCollision& c) { return c.key; } )
 */
template<typename T, typename KeyFn>
static void radixSort(const std::vector<T>& data,
               std::vector<unsigned>& sorted_order,
               unsigned size,
               KeyFn key)
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

        // Count
        for (unsigned idx : sorted_order)
            ++count[(key(data[idx]) >> shift) & MASK];

        // Prefix sums
        size_t sum = 0;
        for (int i = 0; i < BUCKETS; ++i)
        {
            size_t c = count[i];
            count[i] = sum;
            sum += c;
        }

        // Scatter
        for (unsigned idx : sorted_order)
            temp[count[(key(data[idx]) >> shift) & MASK]++] = idx;

        sorted_order.swap(temp);
    }
}

/** Golden section search for minimizing a function F(t) over the interval t_start, t_end.
 * @returns the minimum value t_min
*/
template<typename T, typename F>
static T goldenSectionSearch(const T& t_start, const T& t_end, F&& f, Real tol = 1e-3)
{
    constexpr Real phi = 1.61803398875;
    constexpr Real r = 1 / phi;
    constexpr Real s = 1 - r;

    Real a0 = 0, a1 = s, a2 = r, a3 = 1;

    T t0 = t_start;
    T t1 = lerp(t_start, t_end, a1);
    T t2 = lerp(t_start, t_end, a2);
    T t3 = t_end;

    Real f0 = f(t0), f1 = f(t1), f2 = f(t2), f3 = f(t3);

    while (_gssWidth(t3 - t0) > tol)
    {
        if (std::min(f0, f1) < std::min(f2, f3))
        {
            a3 = a2; t3 = t2; f3 = f2;
            a2 = a1; t2 = t1; f2 = f1;
            a1 = r * a2 + s * a0;
            t1 = lerp(t_start, t_end, a1);
            f1 = f(t1);
        }
        else
        {
            a0 = a1; t0 = t1; f0 = f1;
            a1 = a2; t1 = t2; f1 = f2;
            a2 = r * a1 + s * a3;
            t2 = lerp(t_start, t_end, a2);
            f2 = f(t2);
        }
    }

    T tmid = lerp(t_start, t_end, 0.5 * (a0 + a3));
    Real fmid = f(tmid);

    if (f0 < fmid && f0 < f3) return t0;
    else if (fmid < f3)       return tmid;
    else                      return t3;
}

// helper: "size" of an interval, works for both Real and Vec3r
private:
inline static Real _gssWidth(Real x)        { return std::abs(x); }
inline static Real _gssWidth(const Vec3r& x){ return x.norm(); }


};