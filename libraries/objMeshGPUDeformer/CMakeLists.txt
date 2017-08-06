include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	)

set(OBJMESHBGPUDEFORMER_SRCS
	objMeshGPUDeformer.cpp
	objMeshGPUDeformer_coarseToFine.cpp
	objMeshGPUDeformer_coarseToFine_fbo.cpp
	objMeshGPUDeformer_uUq.cpp
	objMeshGPUDeformer_uUq_fbo.cpp
	objMeshGPUDeformer_uUq_pbuffer.cpp
	objMeshGPUDeformer_uUq_setLighting.cpp
	vbo.cpp
	)

add_library(objMeshGPUDeformer ${OBJMESHBGPUDEFORMER_SRCS})
target_link_libraries(objMeshGPUDeformer objMesh lighting ${CG_LIBRARY} ${CG_GL_LIBRARY})
