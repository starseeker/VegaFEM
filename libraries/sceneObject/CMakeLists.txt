include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(SCOBJ_SRCS
  sceneObject.cpp
  sceneObject6DOF.cpp
  sceneObjectDeformable.cpp
  sceneObjectDeformable6DOF.cpp
  sceneObjectWithRestPosition.cpp
)
if(FOUND_CG)
	set(SCOBJS_SRCS
		${SCOBJS_SRCS}
		sceneObjectDeformableGPU.cpp
		sceneObjectDeformableGPU6DOF.cpp
		)
endif(FOUND_CG)

add_library(sceneObject ${SCOBJ_SRCS})
target_link_libraries(sceneObject objMesh lighting minivector)
