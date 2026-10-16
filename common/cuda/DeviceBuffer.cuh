#pragma once

#include "cuda_runtime.h"

template <typename T>
class DeviceBuffer
{
private:
    T* _ptr = nullptr;
    _sizet _size = 0;

public:
    DeviceBuffer() = default;

    explicit DeviceBuffer(_sizet count)
    {
        allocate(count);
    }

    ~DeviceBuffer()
    {
        reset();
    }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    DeviceBuffer(DeviceBuffer&& other) noexcept
        : _ptr(std::exchange(other._ptr, nullptr))
        , _size(std::exchange(other._size, 0))
    {
    }

    DeviceBuffer& operator=(DeviceBuffer&& other) noexcept
    {
        if (this != &other)
        {
            reset();

            _ptr = std::exchange(other._ptr, nullptr);
            _size = std::exchange(other._size, 0);
        }

        return *this;
    }

    void allocate(_sizet count)
    {
        if (count == _size)
            return;

        reset();

        if (count == 0)
            return;

        T* ptr = nullptr;

        // In production, check the CUDA error.
        cudaMalloc(reinterpret_cast<void**>(&ptr), count * sizeof(T));

        _ptr = ptr;
        _size = count;
    }

    void reset() noexcept
    {
        if (_ptr)
        {
            cudaFree(_ptr);
            _ptr = nullptr;
        }

        _size = 0;
    }

    void upload(const T* src, _sizet count)
    {
        if (count > _size)
            throw std::runtime_error("DeviceBuffer capacity exceeded");

        if (count == 0)
            return;

        cudaMemcpy(_ptr, src, count * sizeof(T), cudaMemcpyHostToDevice);
    }

    void download(T* dst, _sizet count) const
    {
        if (count > _size)
            throw std::runtime_error("DeviceBuffer capacity exceeded");

        if (count == 0)
            return;

        cudaMemcpy(dst, _ptr, count * sizeof(T), cudaMemcpyDeviceToHost);
    }

    T* data() { return _ptr; }
    const T* data() const { return _ptr; }

    _sizet size() const { return _size; }
};