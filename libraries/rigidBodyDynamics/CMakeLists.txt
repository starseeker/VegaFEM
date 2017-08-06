include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}
	)

add_library(rigidBodyDynamics rigidBody.cpp rigidBody_generalTensor.cpp)
target_link_libraries(rigidBodyDynamics quaternion)
