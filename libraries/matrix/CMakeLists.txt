include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(MATRIX_SRCS
	matrixProjection.cpp
	matrixExp.cpp
	matrixBLAS.cpp
	matrix.cpp
	matrixPCA.cpp
  matrixLAPACK.cpp
	)

add_library(matrix ${MATRIX_SRCS})
target_link_libraries(matrix matrixIO ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
