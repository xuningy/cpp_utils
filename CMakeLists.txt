cmake_minimum_required(VERSION 2.8.3)

project(cpp_utils)

find_package(catkin REQUIRED)
find_package(Eigen3 REQUIRED)
catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS
)

catkin_package(
	INCLUDE_DIRS include
	LIBRARIES ${PROJECT_NAME}
)

include_directories(include ${catkin_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
link_directories(${catkin_LIBRARY_DIRS} ${EIGEN3_INCLUDE_DIR})
add_definitions(${catkin_DEFINITIONS})

add_library(${PROJECT_NAME} src/${PROJECT_NAME}/${PROJECT_NAME}.cpp)
target_link_libraries(${PROJECT_NAME}
	${catkin_LIBRARIES}
	stdc++fs
)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

set_target_properties(${PROJECT_NAME}
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
)

install(DIRECTORY include ${PROJECT_NAME}/
DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
FILES_MATCHING PATTERN "*.h"
)
