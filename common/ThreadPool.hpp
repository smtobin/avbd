#pragma once

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>
#include <algorithm>
#include <immintrin.h>   // _mm_pause()

class ThreadPool
{
public:

    explicit ThreadPool(unsigned threads)
        : _numThreads(threads)
    {
        _workers.reserve(_numThreads);

        for (unsigned t = 0; t < _numThreads; ++t)
            _workers.emplace_back([this] { workerLoop(); });
    }

    ~ThreadPool()
    {
        shutdown();
    }

    template<class Func>
    void parallelFor(unsigned count, Func&& func)
    {
        if (count == 0)
            return;

        {
            std::lock_guard lock(_mutex);

            _count = count;
            _next.store(0, std::memory_order_relaxed);
            _finished.store(0, std::memory_order_relaxed);

            _work = &func;

            _callback = [](void* ctx,
                           unsigned begin,
                           unsigned end)
            {
                auto& f = *static_cast<Func*>(ctx);

                f(begin, end);
            };

            ++_generation;
        }

        _workCV.notify_all();

        //
        // Spin briefly.
        // Most colors finish before this expires.
        //
        constexpr int SpinCount = 2000;

        for (int i = 0; i < SpinCount; ++i)
        {
            if (_finished.load(std::memory_order_acquire) == _numThreads)
                return;

            _mm_pause();
        }

        //
        // Longer job.
        // Sleep until last worker finishes.
        //
        std::unique_lock lock(_mutex);

        _doneCV.wait(lock, [&]
        {
            return _finished.load(std::memory_order_acquire) == _numThreads;
        });
    }

private:

    static constexpr unsigned BlockSize = 16;

    using Callback =
        void(*)(void*, unsigned, unsigned);

    //------------------------------------------------------------

    std::vector<std::thread> _workers;

    unsigned _numThreads;

    std::mutex _mutex;

    std::condition_variable _workCV;
    std::condition_variable _doneCV;

    bool _shutdown = false;

    uint64_t _generation = 0;

    unsigned _count = 0;

    Callback _callback = nullptr;

    void* _work = nullptr;

    std::atomic<unsigned> _next{0};

    std::atomic<unsigned> _finished{0};

    //------------------------------------------------------------

    void workerLoop()
    {
        uint64_t myGeneration = 0;

        while (true)
        {
            Callback callback;
            void* work;
            unsigned count;

            {
                std::unique_lock lock(_mutex);

                _workCV.wait(lock, [&]
                {
                    return _shutdown || _generation != myGeneration;
                });

                if (_shutdown)
                    return;

                myGeneration = _generation;

                callback = _callback;
                work = _work;
                count = _count;
            }

            //
            // Dynamic scheduling.
            //
            while (true)
            {
                unsigned begin =
                    _next.fetch_add(BlockSize,
                                    std::memory_order_relaxed);

                if (begin >= count)
                    break;

                unsigned end =
                    std::min(begin + BlockSize,
                             count);

                callback(work, begin, end);
            }

            //
            // Barrier.
            //
            if (_finished.fetch_add(1,
                    std::memory_order_acq_rel)
                + 1 == _numThreads)
            {
                std::lock_guard lock(_mutex);
                _doneCV.notify_one();
            }
        }
    }

    //------------------------------------------------------------

    void shutdown()
    {
        {
            std::lock_guard lock(_mutex);
            _shutdown = true;
        }

        _workCV.notify_all();

        for (auto& t : _workers)
            t.join();
    }
};