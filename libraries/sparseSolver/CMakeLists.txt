include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(SOLVER_SRCS
 linearSolver.cpp
 #ARPACKSolver.cpp
 PardisoSolver.cpp
 SPOOLESSolver.cpp
 SPOOLESSolverMT.cpp
 CGSolver.cpp
 invMKSolver.cpp
)

add_library(sparseSolver ${SOLVER_SRCS})
target_link_libraries(sparseSolver sparseMatrix)
