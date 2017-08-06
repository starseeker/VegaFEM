#FindOpenGL
find_package(OpenGL)
# #FindGLUT
find_package(GLUT)
# #FindGLEW
find_package(GLEW)
# #FindCg
find_package(Cg)
#FindLAPACK
find_package(LAPACK)
#FindBLAS
find_package(BLAS)
#FindThreads
find_package(Threads REQUIRED)

include_directories(
	${OPENGL_INCLUDE_DIR}
	${GLUT_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/libraries/include
	)


set(VEGA_CORE_LIBRARIES
	clothBW
	configFile
	corotationalLinearFEM
	elasticForceModel
	forceModel
	getopts
	graph
	hashTable
	imageIO
	insertRows
	integrator
	integratorSparse
	isotropicHyperelasticFEM
	loadList
	massSpringSystem
	matrixIO
	minivector
	objMesh
	performanceCounter
	polarDecomposition
	quaternion
	rigidBodyDynamics
	sparseMatrix
	sparseSolver
	stvk
	volumetricMesh
	)

if(OPENGL_FOUND AND FOUND_CG)
	list(APPEND VEGA_CORE_LIBRARIES objMeshGPUDeformer)
endif(OPENGL_FOUND AND FOUND_CG)

if(OPENGL_FOUND AND GLUT_FOUND)
	list(APPEND VEGA_CORE_LIBRARIES
		glslPhong
		openGLHelper
		camera
		lighting
		sceneObject
		renderVolumetricMesh
		)
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

