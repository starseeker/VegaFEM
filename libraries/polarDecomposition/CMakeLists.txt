include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(POLAR_DECOMPOSITION_SRCS
  polarDecomposition.cpp
  polarDecompositionGradient.cpp
)

add_library(polarDecomposition ${POLAR_DECOMPOSITION_SRCS})
