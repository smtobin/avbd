#pragma once

#include <atomic>

class SpinBarrier
{
private:
    alignas(64) std::atomic<unsigned> _arrived = 0;
    alignas(64) std::atomic<unsigned> _generation = 0;  
    
    unsigned _thread_count;

public:
    SpinBarrier()
        : _thread_count(1)
    {}

    SpinBarrier(unsigned thread_count)
        : _thread_count(thread_count)
    {}

    void arrive_and_wait()
    {
        unsigned gen = _generation.load(std::memory_order_acquire);

        // increment arrived
        // if this is same as the number of threads, increment the generation and return
        if (_arrived.fetch_add(1) == _thread_count - 1)
        {
            _arrived.store(0, std::memory_order_relaxed);
            _generation.fetch_add(1, std::memory_order_release);
            return;
        }

        // spin until generation is increased
        while (_generation.load(std::memory_order_acquire) == gen)
            _mm_pause();
    }
};