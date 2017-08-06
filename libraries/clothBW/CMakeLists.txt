include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(CLOTHBW_SRCS
	clothBW.cpp
	clothBWFromObjMesh.cpp
	clothBWMT.cpp
)

add_library(clothBW ${CLOTHBW_SRCS})
target_link_libraries(clothBW minivector objMesh sparseMatrix)
