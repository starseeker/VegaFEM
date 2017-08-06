include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(SM_SRCS
  sparseMatrix.cpp
  sparseMatrixMT.cpp
)

add_library(sparseMatrix ${SM_SRCS})
