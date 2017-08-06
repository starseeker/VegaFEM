include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
	)

set(integratorsparse_srcs
	centralDifferencesSparse.cpp
	eulerSparse.cpp
	implicitBackwardEulerSparse.cpp
	implicitNewmarkSparse.cpp
	integratorBaseSparse.cpp
	)

add_library(integratorSparse ${integratorsparse_srcs})
target_link_libraries(integratorSparse integrator performanceCounter insertRows sparseSolver forceModel)
