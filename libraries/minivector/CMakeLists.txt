include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(MINIVEC_SRCS
   mat3d.cpp
   vec3d.cpp
   vec2d.cpp
   eig3.cpp
)

add_library(minivector ${MINIVEC_SRCS})
