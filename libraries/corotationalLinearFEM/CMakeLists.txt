include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(corotationalLinearFEM corotationalLinearFEM.cpp corotationalLinearFEMMT.cpp)
target_link_libraries(corotationalLinearFEM polarDecomposition volumetricMesh sparseMatrix)
