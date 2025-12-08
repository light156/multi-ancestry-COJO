#pragma once

#if defined(_OPENMP)
    #include <omp.h>
    #define HAS_OPENMP 1
    #define OMP_PARALLEL_FOR _Pragma("omp parallel for schedule(dynamic)")
#else
    #define HAS_OPENMP 0
    #define OMP_PARALLEL_FOR
    inline void omp_set_num_threads(int) {}
#endif
