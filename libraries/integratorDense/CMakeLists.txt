
set(INTEGRATORDENSE_SRCS
	centralDifferencesDense.cpp
	implicitBackwardEulerDense.cpp
	implicitNewmarkDense.cpp
	integratorBaseDense.cpp
	)

add_library(integratordense ${INTEGRATORDENSE_SRCS})
target_link_libraries(integratordense matrix integrator performanceCounter reducedForceModel)
