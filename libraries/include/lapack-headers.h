#if defined(WIN32) || defined(_WIN32) || defined(linux)
  #include "mkl_cblas.h"
  #include "mkl_types.h"
  #include "mkl_lapack.h"
  #include "mkl_blas.h"
#elif defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
  #if MAC_OS_X_VERSION_MIN_REQUIRED == MAC_OS_X_VERSION_10_9
  // On one of our MacBook with Mac OS X 10.9, the following headers are required to find CBLAS routines.
    #include <vecLib/cblas.h>
    #include <vecLib/clapack.h>
  #endif
#endif

