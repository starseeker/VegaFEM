include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(STVK_SRCS
 StVKCubeABCD.cpp
 StVKElementABCD.cpp
 StVKElementABCDLoader.cpp
 StVKHessianTensor.cpp
 StVKInternalForces.cpp
 StVKInternalForcesMT.cpp
 StVKStiffnessMatrix.cpp
 StVKStiffnessMatrixMT.cpp
 StVKTetABCD.cpp
 StVKTetHighMemoryABCD.cpp
)

add_library(stvk ${STVK_SRCS})
target_link_libraries(stvk minivector volumetricMesh sparseMatrix)
