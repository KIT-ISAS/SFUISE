cmake_minimum_required(VERSION 2.8.3)
project(cf_msgs)

find_package(catkin REQUIRED message_generation std_msgs)


add_message_files(
  FILES
  Gyro.msg
  Accel.msg
  Tdoa.msg
  Twr.msg
  Tof.msg
  Baro.msg
  Flow.msg
)

generate_messages(
   DEPENDENCIES
   std_msgs  # Or other packages containing msgs
)

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES cf_msgs
  CATKIN_DEPENDS std_msgs message_runtime
#  DEPENDS system_lib
)

include_directories(
# include
# ${catkin_INCLUDE_DIRS}
)
