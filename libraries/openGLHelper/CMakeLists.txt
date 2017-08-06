include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(OPENGLHELPER_SRCS
	openGLHelper.cpp
	printBitmap.cpp
	saveScreenShot.cpp
	)

add_library(openGLHelper ${OPENGLHELPER_SRCS})
target_link_libraries(openGLHelper matrixIO imageIO)
