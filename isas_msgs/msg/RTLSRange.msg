uint16 id # anchor id
geometry_msgs/Point pos # anchor position
float32 pr # position rate for the quality of the estimated anchor position
float32 range # current distance to the anchor in m
float32 mean # mean value of distances since reset
float32 var # variance of the distance
float32 fpp # first path power
float32 rxp # received power
float32 csn # channel standard noise
float32 cmn # channel maximum noise
uint8 toc # time out counter, increases if the range measurement fails
uint8 ra # 1 denotes valid data and has been used for the estimation of position, 0 denotes old data, the measurement didn't success
