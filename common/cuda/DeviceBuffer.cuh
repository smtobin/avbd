#pragma once

#include "cuda_runtime.h"

#include "common/common.hpp"
#include "common/cuda/Vec3.cuh"
#include "common/cuda/Quaternion.cuh"

#include <vector>

template <typename T>
class DeviceBuffer
{
private:
    T* _ptr = nullptr;
    size_t _size = 0;

public:
    DeviceBuffer() = default;

    explicit DeviceBuffer(size_t count)
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

    void allocate(size_t count)
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

    void upload(const T* src, size_t count)
    {
        if (count > _size)
            throw std::runtime_error("DeviceBuffer capacity exceeded");

        if (count == 0)
            return;

        cudaMemcpy(_ptr, src, count * sizeof(T), cudaMemcpyHostToDevice);
    }

    void download(T* dst, size_t count) const
    {
        if (count > _size)
            throw std::runtime_error("DeviceBuffer capacity exceeded");

        if (count == 0)
            return;

        cudaMemcpy(dst, _ptr, count * sizeof(T), cudaMemcpyDeviceToHost);
    }

    T* data() { return _ptr; }
    const T* data() const { return _ptr; }

    size_t size() const { return _size; }
};

// add specific overloads for uploading Vec3r and Quaternion
void DeviceBuffer<DVec3>::upload(const Vec3r* src, size_t count)
{
    if (count > _size)
            throw std::runtime_error("DeviceBuffer capacity exceeded");

    if (count == 0)
        return;

    // safe approach - convert Eigen Vec3r to device DVec3, then copy
    std::vector<DVec3> staging(count);

    for (size_t i = 0; i < count; ++i)
    {
        staging[i] = DVec3(src[i].x(), src[i].y(), src[i].z());
    }

    const cudaError_t err = cudaMemcpy(
        _ptr,
        staging.data(),
        count * sizeof(DVec3),
        cudaMemcpyHostToDevice
    );

    if (err != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(err));
}

void DeviceBuffer<DQuaternion>::upload(const Quaternion* src, size_t count)
{
    if (count > _size)
            throw std::runtime_error("DeviceBuffer capacity exceeded");

    if (count == 0)
        return;

    // safe approach - convert Eigen Quaternion to device DQuaternion, then copy
    std::vector<DQuaternion> staging(count);

    for (size_t i = 0; i < count; ++i)
    {
        staging[i] = DQuaternion(src[i].w(), src[i].x(), src[i].y(), src[i].z());
    }

    const cudaError_t err = cudaMemcpy(
        _ptr,
        staging.data(),
        count * sizeof(DQuaternion),
        cudaMemcpyHostToDevice
    );

    if (err != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(err));
}