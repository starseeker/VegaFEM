include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(glslPhong glslPhong.cpp)

target_link_libraries(glslPhong ${GLEW_LIBRARIES})
