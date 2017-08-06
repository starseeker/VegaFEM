include_directories(
   ${CMAKE_CURRENT_SOURCE_DIR}
)

set(IHFEM_SRCS
 isotropicMaterial.cpp
 MooneyRivlinIsotropicMaterial.cpp
 neoHookeanIsotropicMaterial.cpp
 StVKIsotropicMaterial.cpp
 homogeneousMooneyRivlinIsotropicMaterial.cpp
 homogeneousStVKIsotropicMaterial.cpp
 homogeneousNeoHookeanIsotropicMaterial.cpp
 isotropicHyperelasticFEM.cpp
 isotropicHyperelasticFEMMT.cpp
 isotropicMaterialWithCompressionResistance.cpp
)

add_library(isotropicHyperelasticFEM ${IHFEM_SRCS})
target_link_libraries(isotropicHyperelasticFEM minivector volumetricMesh sparseMatrix)
