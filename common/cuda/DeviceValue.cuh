
#pragma once

#include <cuda_runtime.h>

#include <type_traits>
#include <stdexcept>
#include <string>
#include <utility>

template <typename T>
class DeviceValue {
    static_assert(std::is_trivially_copyable<T>::value, "DeviceValue<T> requires T to be trivially copyable.");
    static_assert(std::is_default_constructible<T>::value, "DeviceValue<T> requires T to be default constructible.");
private:
    T* _ptr = nullptr;

public:
    using value_type = T;

    // Empty, unallocated value.
    DeviceValue() noexcept = default;

    // Allocate storage and initialize it from a host value.
    explicit DeviceValue(const T& initial_value) 
    {
        allocate();

        try 
        {
            set(initial_value);
        } catch (...) 
        {
            _releaseNoThrow();
            throw;
        }
    }

    ~DeviceValue() noexcept { _releaseNoThrow(); }

    // Device memory has unique ownership.
    DeviceValue(const DeviceValue&) = delete;
    DeviceValue& operator=(const DeviceValue&) = delete;

    // Transfer ownership without copying device contents.
    DeviceValue(DeviceValue&& other) noexcept
        : _ptr(std::exchange(other._ptr, nullptr)) {}

    DeviceValue& operator=(DeviceValue&& other) noexcept {
        if (this != &other) 
        {
            _releaseNoThrow();
            _ptr = std::exchange(other._ptr, nullptr);
        }
        return *this;
    }

    // Allocate one T. The value is uninitialized until written.
    void allocate() 
    {
        if (_ptr != nullptr)
            return;

        T* ptr = nullptr;
        checkCuda(cudaMalloc(reinterpret_cast<void**>(&ptr), sizeof(T)), "cudaMalloc");
        _ptr = ptr;
    }

    // Free the allocation. Throws if cudaFree reports an error.
    void reset() 
    {
        if (_ptr == nullptr)
            return;

        checkCuda(cudaFree(_ptr), "cudaFree");
        _ptr = nullptr;
    }

    bool allocated() const noexcept 
    {
        return _ptr != nullptr;
    }

    // Pointer to the device-resident value.
    // Pass this pointer to kernels; do not dereference it on the CPU.
    T* data() noexcept { return _ptr; }

    const T* data() const noexcept { return _ptr; }

    // Synchronous host -> device transfer.
    void set(const T& value) 
    {
        _requireAllocated();

        checkCuda(
            cudaMemcpy(
                _ptr,
                &value,
                sizeof(T),
                cudaMemcpyHostToDevice),
            "cudaMemcpy (host to device)");
    }

    // Synchronous device -> host transfer.
    T get() const 
    {
        _requireAllocated();

        T result{};
        checkCuda(
            cudaMemcpy(
                &result,
                _ptr,
                sizeof(T),
                cudaMemcpyDeviceToHost),
            "cudaMemcpy (device to host)");

        return result;
    }

    // Asynchronous host -> device transfer.
    //
    // The source pointer must remain valid until the transfer completes.
    // Pinned host memory is recommended for genuinely asynchronous copies.
    void setAsync(const T* host_source, cudaStream_t stream = nullptr) 
    {
        _requireAllocated();

        if (host_source == nullptr)
            throw std::invalid_argument(
                "DeviceValue::setAsync received a null source.");

        checkCuda(
            cudaMemcpyAsync(
                _ptr,
                host_source,
                sizeof(T),
                cudaMemcpyHostToDevice,
                stream),
            "cudaMemcpyAsync (host to device)");
    }

    // Asynchronous device -> host transfer.
    //
    // The destination pointer must remain valid until the transfer
    // completes. Synchronize the stream before reading the destination.
    void getAsync(T* host_destination, cudaStream_t stream = nullptr) const 
    {
        _requireAllocated();

        if (host_destination == nullptr)
            throw std::invalid_argument(
                "DeviceValue::getAsync received a null destination.");

        checkCuda(
            cudaMemcpyAsync(
                host_destination,
                _ptr,
                sizeof(T),
                cudaMemcpyDeviceToHost,
                stream),
            "cudaMemcpyAsync (device to host)");
    }

private:

    void _requireAllocated() const 
    {
        if (_ptr == nullptr)
            throw std::logic_error(
                "DeviceValue operation requires allocated storage.");
    }

    static void checkCuda(cudaError_t error, const char* operation) 
    {
        if (error != cudaSuccess) {
            throw std::runtime_error(
                std::string(operation) + " failed: " +
                cudaGetErrorString(error));
        }
    }

    // Destructors and move assignment must not throw.
    void _releaseNoThrow() noexcept 
    {
        if (_ptr != nullptr) {
            (void)cudaFree(_ptr);
            _ptr = nullptr;
        }
    }
};
