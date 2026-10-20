#pragma once

#ifndef HD
    #ifdef __CUDACC__
        #define HD __host__ __device__
    #else
        #define HD
    #endif
#endif

template <typename Real>
HD Real absValue(Real x)
{
    return x < Real(0) ? -x : x;
}