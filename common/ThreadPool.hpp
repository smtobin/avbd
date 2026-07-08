#pragma once

#include <vector>
#include <functional>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>

class ThreadPool {
public:
    explicit ThreadPool(unsigned num_threads) {
        // create workers that will wait for work
        for (unsigned t = 0; t < num_threads; t++)
            workers.emplace_back([this, t] { workerLoop(t); });
    }

    // Dispatch a parallel-for over [0, n) and block until complete
    template <typename Func>
    void parallelFor(unsigned n, Func&& f) {
        if (n == 0) return;

        work = [&f, n, num_threads = workers.size()](unsigned t) {
            // divvy up the work equally between the worker threads
            unsigned chunk = (n + num_threads - 1) / num_threads;
            unsigned start = t * chunk;
            unsigned end   = std::min(start + chunk, n);

            for (unsigned i = start; i < end; i++) 
                f(i);
        };
        dispatch();
    }

    ~ThreadPool() { shutdown(); }

private:
    std::vector<std::thread> workers;
    std::function<void(unsigned)> work;
    std::atomic<unsigned> done_count{0};
    std::atomic<bool> stop{false};
    std::atomic<long unsigned> generation{0};
    std::condition_variable cv;
    std::mutex mtx;

    void dispatch() {
        // increment generation to signal to other threads that new work is ready
        ++generation;

        // keep track of how many threads are done
        done_count = 0;
        { std::lock_guard<std::mutex> lk(mtx); cv.notify_all(); }
        while (done_count.load(std::memory_order_acquire) < workers.size()) {
            std::this_thread::yield();  // spin-wait — see note below
        }
    }

    void workerLoop(unsigned t) {
        // track the number of times the CV has been notified
        // used to unlock and start doing work
        long unsigned last_gen = 0;

        while (!stop.load(std::memory_order_acquire)) {
            // wait for work
            std::unique_lock<std::mutex> lk(mtx);
            cv.wait(lk, [&] { return generation.load() != last_gen || stop.load(); });
            lk.unlock();

            // quit if stop is true at any point
            if (stop.load()) 
                return;
            
            
            last_gen = generation.load();
            
            // perform the work
            work(t);

            // increment done count when done
            done_count.fetch_add(1, std::memory_order_release);
        }
    }

    void shutdown() {
        stop = true;
        { std::lock_guard<std::mutex> lk(mtx); cv.notify_all(); }
        for (auto& w : workers) w.join();
    }
};