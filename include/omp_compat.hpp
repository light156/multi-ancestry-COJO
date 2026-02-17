#pragma once

#if defined(_OPENMP)
    #include <omp.h>
    #define HAS_OPENMP 1
#else
    #define HAS_OPENMP 0
    inline void omp_set_num_threads(int) {}
    inline int omp_get_thread_num() { return 0; }
#endif
