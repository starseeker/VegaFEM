#ifdef WIN32
#ifdef _DEBUG
#pragma comment(lib, "glew32d.lib")
#pragma comment(lib, "freeglutd.lib")
// image libraries
// #pragma comment(lib, "tiffd.lib")
// #pragma comment(lib, "jpeg.lib")
// #pragma comment(lib, "libpng16.lib")
#else
#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "freeglut.lib")
// image libraries
// #pragma comment(lib, "tiff.lib")
// #pragma comment(lib, "jpeg.lib")
// #pragma comment(lib, "libpng16.lib")
#endif
#endif

// mkl
// #pragma comment(lib, "mkl_core.lib")
// #pragma comment(lib, "mkl_intel_thread.lib")
// #pragma comment(lib, "mkl_intel_lp64.lib")
// #pragma comment(lib, "libiomp5md.lib")
// Cg
// #pragma comment(lib, "cg.lib")
// #pragma comment(lib, "cgGL.lib")
