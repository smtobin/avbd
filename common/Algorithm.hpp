#pragma once

#include "common/common.hpp"
#include "common/WorkerThreadContext.hpp"

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

/** Parallel radix sort on uint64 keys.
 * @param w_ctx this worker thread context
 * @param all_contexts all the worker contexts
 * @param combined_offsets shared scratch memory with size = num_threads * 256
 * @param unsorted the list of unsorted keys
 * @param sorted_order (output) the lindices of the unsorted array, arranged in the sorted order
 * @param temp shared scratch buffer, same size as sorted_order
 * @param size the size of the vector (the unsorted vector may have extra padding)
 */
static void radixSort_Parallel(
    WorkerThreadContext& w_ctx,
    const std::vector<WorkerThreadContext>& all_contexts,
    std::vector<unsigned>& combined_offsets, 
    const std::vector<uint64_t>& unsorted,
    std::vector<unsigned>& sorted_order,
    std::vector<unsigned>& temp,
    unsigned size)
{
    constexpr int BITS = 8;
    constexpr int BUCKETS = 1 << BITS;
    constexpr int MASK = BUCKETS - 1;
    const unsigned num_threads = WorkerThreadContext::NUM_THREADS;

    // allocated space
    if (w_ctx.idx == 0)
    {
        sorted_order.resize(size);
        std::iota(sorted_order.begin(), sorted_order.end(), 0);
        temp.resize(size);
        combined_offsets.resize(num_threads * BUCKETS);
    }
    w_ctx.barrier->arrive_and_wait();

    auto [start, end] = w_ctx.computeStartEnd(size);

    for (int shift = 0; shift < 64; shift += BITS)
    {
        // local histogram
        auto& local_count = w_ctx.RadixSortContext.local_count;
        std::fill(std::begin(local_count), std::end(local_count), 0);

        for (unsigned i = start; i < end; i++)
        {
            ++local_count[(unsorted[sorted_order[i]] >> shift) & MASK];
        }

        w_ctx.barrier->arrive_and_wait();
    

        // leader builds (bucket, thread) ==> global offset
        if (w_ctx.idx == 0)
        {
            unsigned running = 0;
            for (int b = 0; b < BUCKETS; b++)
            {
                for (unsigned t = 0; t < num_threads; t++)
                {
                    combined_offsets[t*BUCKETS + b] = running;
                    running += all_contexts[t].RadixSortContext.local_count[b];
                }
            }
        }

        w_ctx.barrier->arrive_and_wait();

        // parallel scatter
        unsigned cursor[BUCKETS];
        // copy from shared combined_offsets into local cursor
        std::copy(
            &combined_offsets[w_ctx.idx * BUCKETS],
            &combined_offsets[w_ctx.idx * BUCKETS] + BUCKETS,
            cursor);
        for (unsigned i = start; i < end; i++)
        {
            unsigned idx = sorted_order[i];
            unsigned bucket = (unsorted[idx] >> shift) & MASK;
            temp[cursor[bucket]++] = idx;
        }

        w_ctx.barrier->arrive_and_wait();

        // swap buffers (with single thread)
        if (w_ctx.idx == 0)
            sorted_order.swap(temp);
        
        w_ctx.barrier->arrive_and_wait();
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

/** Parallel prefix sum
 * @param w_ctx the worker context
 * @param counts the counts for each item in the list
 * @param out the output offsets for each item in the list
 * @param n the number of items in the list
 * @param chunk_totals shared scratch memory - size = num threads
 */
static void parallelPrefixSum(
    WorkerThreadContext& w_ctx,
    const unsigned* counts,
    unsigned* out,
    unsigned n,
    std::vector<unsigned>& chunk_totals
)
{
    auto [start, end] = w_ctx.computeStartEnd(n);

    // each thread goes through its range and accumulates thread-local offsets
    unsigned running = 0;
    for (unsigned i = start; i < end; i++)
    {
        out[i] = running;
        running += counts[i];
    }
    // store the total in chunk_totals for each thread
    chunk_totals[w_ctx.idx] = running;
    w_ctx.barrier->arrive_and_wait();

    // main thread will accumulate the chunk_totals
    // chunk_totals will now store the base offset for each thread
    if (w_ctx.idx == 0)
    {
        unsigned base = 0;
        for (unsigned t = 0; t < WorkerThreadContext::NUM_THREADS; t++)
        {
            unsigned tmp = chunk_totals[t];
            chunk_totals[t] = base;
            base += tmp;
        }
        out[n] = base; // grand total
    }
    w_ctx.barrier->arrive_and_wait();

    unsigned base = chunk_totals[w_ctx.idx];
    for (unsigned i = start; i < end; i++)
    {
        // increment the offsets based on the global start
        out[i] += base;
    }
    w_ctx.barrier->arrive_and_wait();
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