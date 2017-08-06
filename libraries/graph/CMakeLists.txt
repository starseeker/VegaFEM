include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(graph graph.cpp)
target_link_libraries(graph sparseMatrix)
