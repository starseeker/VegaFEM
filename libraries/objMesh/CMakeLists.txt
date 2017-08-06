include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(OBJ_SRCS
  objMesh.cpp
  objMesh-disjointSet.cpp
  objMeshBinaryLoader.cpp
  objMeshEncode.cpp
  objMeshGraph.cpp
  objMeshOctree.cpp
  objMeshOffsetVoxels.cpp
  objMeshOrientable.cpp
  simpleSphere.cpp
  tribox3.cpp
)

if(OPENGL_FOUND AND GLUT_FOUND)
 list(APPEND OBJ_SRCS
  boundingBox.cpp
  objMeshRender.cpp
  octree.cpp
  triangle.cpp
  )
endif()

add_library(objMesh ${OBJ_SRCS})
target_link_libraries(objMesh minivector matrixIO imageIO sparseMatrix graph)

if(OPENGL_FOUND AND GLUT_FOUND)
  target_link_libraries(objMesh ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
endif()
