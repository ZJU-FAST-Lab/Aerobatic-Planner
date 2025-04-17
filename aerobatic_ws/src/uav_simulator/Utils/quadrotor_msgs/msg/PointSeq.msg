Header header

# the action command for trajectory server.
uint8 POINT_TYPE_0                =   0
uint8 POINT_TYPE_1                =   1
uint8 POINT_TYPE_2                =   2

uint32 length

# the polynomial coecfficients of the trajectory.
float64[] pos
float64[] vel
float64[] acc
uint8[] type