#include "vbo.h"
// #include <iostream>
// using namespace std;
#if defined(_WIN32) || defined(WIN32)

#include "glh_extensions.h"

//extern PFNGLBINDBUFFERARBPROC glBindBufferARB;
//extern PFNGLBUFFERDATAARBPROC glBufferDataARB;
//extern PFNGLGENBUFFERSARBPROC glGenBuffersARB;
//extern PFNGLDELETEBUFFERSARBPROC glDeleteBuffersARB;

bool InitializeVBOs(void)
{
  // cout << "Initializing VBOs (in Windows)." << endl;
  if (glBindBufferARB == NULL)
    glBindBufferARB = (PFNGLBINDBUFFERARBPROC) wglGetProcAddress("glBindBufferARB");

  if (glBufferDataARB == NULL)
    glBufferDataARB = (PFNGLBUFFERDATAARBPROC) wglGetProcAddress("glBufferDataARB");

  if (glGenBuffersARB == NULL)
    glGenBuffersARB = (PFNGLGENBUFFERSARBPROC) wglGetProcAddress("glGenBuffersARB");

  if (glDeleteBuffersARB == NULL)
    glDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC) wglGetProcAddress("glDeleteBuffersARB");

  // cout << "glGenBuffersARB: " << glGenBuffersARB << endl;
  // cout << "return " << (glBindBufferARB && glBufferDataARB && glGenBuffersARB && glDeleteBuffersARB) << endl;
  return (glBindBufferARB && glBufferDataARB && glGenBuffersARB && glDeleteBuffersARB);
}

#else
  bool InitializeVBOs(void) { return true; }
#endif

