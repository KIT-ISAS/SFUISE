cmake_minimum_required(VERSION 3.0.2)
project(isas_msgs)

find_package(catkin REQUIRED message_generation std_msgs geometry_msgs)

add_message_files(
  FILES
  RTLSRange.msg
  RTLSStick.msg
  AnchorPosition.msg
  Anchorlist.msg
)

generate_messages(
 DEPENDENCIES
 std_msgs
 geometry_msgs
)

catkin_package(
  CATKIN_DEPENDS std_msgs message_runtime
)

