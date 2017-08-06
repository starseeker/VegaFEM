include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	)

add_library(reducedForceModel reducedForceModel.cpp reducedForceModelWithHessian.cpp)
