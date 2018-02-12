cmake_minimum_required(VERSION 3.1)

#FindOpenGL
find_package(OpenGL QUIET)
# #FindGLUT
find_package(GLUT QUIET)
# #FindGLEW
find_package(GLEW QUIET)
# #FindCg
find_package(Cg QUIET)
#FindBLAS
find_package(BLAS QUIET)
if(BLAS_FOUND)
	#FindLAPACK
	find_package(LAPACK QUIET)
endif(BLAS_FOUND)
#FindThreads
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads QUIET)

include_directories(
	${OPENGL_INCLUDE_DIR}
	${GLUT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/libraries/include
	)

set(VEGA_CORE_LIBRARIES
	configFile
	forceModel
	getopts
	graph
	hashTable
	imageIO
	insertRows
	integrator
	integratorSparse
	loadList
	matrixIO
	minivector
	objMesh
	performanceCounter
	polarDecomposition
	quaternion
	rigidBodyDynamics
	sparseMatrix
	sparseSolver
	volumetricMesh
	)

if(CMAKE_USE_PTHREADS_INIT)
	list(APPEND VEGA_CORE_LIBRARIES
		clothBW
		corotationalLinearFEM
		elasticForceModel
		isotropicHyperelasticFEM
		massSpringSystem
		stvk
		)
endif(CMAKE_USE_PTHREADS_INIT)

if(OPENGL_FOUND AND FOUND_CG)
	list(APPEND VEGA_CORE_LIBRARIES objMeshGPUDeformer)
endif(OPENGL_FOUND AND FOUND_CG)

if(OPENGL_FOUND AND GLUT_FOUND)
	list(APPEND VEGA_CORE_LIBRARIES
		openGLHelper
		camera
		lighting
		sceneObject
		renderVolumetricMesh
		)
	if(GLEW_FOUND)
		list(APPEND VEGA_CORE_LIBRARIES glslPhong)
	endif(GLEW_FOUND)
endif(OPENGL_FOUND AND GLUT_FOUND)

if(LAPACK_FOUND AND BLAS_FOUND)
	list(APPEND VEGA_CORE_LIBRARIES
		reducedElasticForceModel
		reducedForceModel
		reducedStvk
		sceneObjectReduced
		integratorDense
		modalMatrix
		matrix
		)
endif(LAPACK_FOUND AND BLAS_FOUND)

foreach(lib ${VEGA_CORE_LIBRARIES})
	include_directories(libraries/${lib})
endforeach(lib ${VEGA_CORE_LIBRARIES})

foreach(lib ${VEGA_CORE_LIBRARIES})
	add_subdirectory(libraries/${lib})
endforeach(lib ${VEGA_CORE_LIBRARIES})

if(BUILD_UTILITIES)
	add_subdirectory(utilities)
endif(BUILD_UTILITIES)

add_library(vega "${CMAKE_CURRENT_SOURCE_DIR}/libraries/dummy.cpp")
target_link_libraries(vega ${VEGA_CORE_LIBRARIES})
target_include_directories(vega PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libraries)
target_include_directories(vega PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libraries/include)
foreach(lib ${VEGA_CORE_LIBRARIES})
	target_include_directories(vega PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib})
endforeach(lib ${VEGA_CORE_LIBRARIES})