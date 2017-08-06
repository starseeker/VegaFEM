include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	)

set(rstvk_srcs
	StVKReducedInternalForces.cpp
	StVKReducedInternalForcesMT.cpp
	StVKReducedStiffnessMatrix.cpp
	StVKReducedHessianTensor.cpp
	)

add_library(reducedStvk ${rstvk_srcs})
target_link_libraries(reducedStvk stvk modalMatrix matrix)
