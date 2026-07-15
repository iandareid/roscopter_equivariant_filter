cmake_minimum_required(VERSION 3.8)
project(roscopter_eqf)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(rclcpp REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(roscopter REQUIRED)
find_package(roscopter_msgs REQUIRED)
find_package(rosflight_msgs REQUIRED)

include_directories(
  include
  ${EIGEN3_INCLUDE_DIRS}
)

ament_export_include_directories(include ${EIGEN3_INCLUDE_DIRS})
ament_export_dependencies(
  Eigen3
  rclcpp
  geometry_msgs
  sensor_msgs
  roscopter
  roscopter_msgs
  rosflight_msgs
)

add_library(param_manager
  src/param_manager/param_manager.cpp
)
ament_target_dependencies(param_manager rclcpp)
target_include_directories(param_manager PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
ament_export_targets(param_manager HAS_LIBRARY_TARGET)

add_executable(eqf
  src/eqf/roscopter_eqf.cpp
)
ament_target_dependencies(eqf
  Eigen3
  geometry_msgs
  rclcpp
  roscopter
  roscopter_msgs
  rosflight_msgs
  sensor_msgs
)
target_link_libraries(eqf
  roscopter::estimation_lib
  param_manager
)
target_include_directories(eqf PUBLIC include)

install(DIRECTORY include/
  DESTINATION include
)

install(TARGETS
  param_manager
  EXPORT param_manager
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

install(TARGETS
  eqf
  DESTINATION lib/${PROJECT_NAME}
)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  set(ament_cmake_copyright_FOUND TRUE)
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
