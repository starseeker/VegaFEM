include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(INTEGRATOR_SRCS
 integratorBase.cpp
 getIntegratorSolver.cpp
)

add_library(integrator ${INTEGRATOR_SRCS})
target_link_libraries(integrator performanceCounter insertRows sparseSolver forceModel)
