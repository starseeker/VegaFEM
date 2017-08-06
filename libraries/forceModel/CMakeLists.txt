include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(forceModel forceModel.cpp)
target_link_libraries(forceModel sparseMatrix)
