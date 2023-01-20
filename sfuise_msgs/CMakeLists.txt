cmake_minimum_required(VERSION 3.0.2)
project(sfuise_msgs)

find_package(catkin REQUIRED message_generation std_msgs geometry_msgs)

add_message_files(
	FILES
	Spline.msg
	Knot.msg
	Calib.msg
	Offset.msg
	Estimate.msg
)

generate_messages(
 DEPENDENCIES
 geometry_msgs
)

catkin_package(
  CATKIN_DEPENDS message_runtime
)

