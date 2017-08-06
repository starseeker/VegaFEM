include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(EFM_SRCS
 corotationalLinearFEMForceModel.cpp
 massSpringSystemForceModel.cpp
 StVKForceModel.cpp
 isotropicHyperelasticFEMForceModel.cpp
 linearFEMForceModel.cpp
)

add_library(elasticForceModel ${EFM_SRCS})
target_link_libraries(elasticForceModel forceModel corotationalLinearFEM massSpringSystem stvk isotropicHyperelasticFEM)
