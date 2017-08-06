include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

add_library(lighting lighting.cpp)
target_link_libraries(lighting configFile)
