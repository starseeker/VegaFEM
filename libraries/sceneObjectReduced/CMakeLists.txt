include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	)

set(SCENE_OBJREDUCED_SRCS
	sceneObjectReduced.cpp
	sceneObjectReduced6DOF.cpp
	sceneObjectReducedCPU.cpp
	sceneObjectReducedCPU6DOF.cpp
	)
if(CG_FOUND)
	set(SCENEOBJREDUCED_SRCS
		${SCENEOBJREDUCED_SRCS}
		sceneObjectReducedGPU.cpp
		sceneObjectReducedGPU6DOF.cpp
		)
endif(CG_FOUND)

add_library(sceneObjectReduced ${SCENE_OBJREDUCED_SRCS})
target_link_libraries(sceneObjectReduced objMesh lighting minivector modalMatrix sceneObject)
if(CG_FOUND)
	target_link_libraries(sceneObjectReduced objMeshGPUDeformer)
endif(CG_FOUND)

