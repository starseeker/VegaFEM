include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	)

add_library(renderVolumetricMesh renderVolumetricMesh.cpp)
target_link_libraries(renderVolumetricMesh volumetricMesh openGLHelper)
