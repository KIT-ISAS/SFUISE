Header header
uint16 id
uint32 t # timestamp in ms since the modul is turned on
geometry_msgs/Point pos # current position
geometry_msgs/Point mean # mean position
geometry_msgs/Vector3 hrp # heading, roll, pitch
uint8 noga # number of known anchors
uint8 nora # number of anchors used for the current estimation of position
RTLSRange[] ranges # range measurements
