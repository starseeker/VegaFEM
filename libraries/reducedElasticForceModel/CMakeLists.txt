include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	)

set(REDUCED_ELASTICFORCEMODEL_SRCS
	reducedLinearForceModel.cpp
	reducedLinearStVKForceModel.cpp
	reducedMassSpringSystemForceModel.cpp
	reducedMassSpringSystemForceModelWithHessian.cpp
	reducedSpringForceModel.cpp
	reducedStVKForceModel.cpp
	reducedStVKForceModelWithHessian.cpp
	reducedSubspaceStVKForceModel.cpp
	)

add_library(reducedElasticForceModel ${REDUCED_ELASTICFORCEMODEL_SRCS})
target_link_libraries(reducedElasticForceModel reducedForceModel reducedStvk stvk massSpringSystem)

