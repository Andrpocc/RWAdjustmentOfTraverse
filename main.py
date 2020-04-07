from parametric_adjust import parametric_adjustment, adjustment_traverse, to_degrees
from random import normalvariate
from math import sqrt

mb = 9
ms = 2 + 2 * s

input_angles = [to_degrees(206, 14, 53), to_degrees(146, 2, 57), to_degrees(285, 26, 59),
               to_degrees(162, 44, 31), to_degrees(96, 6, 52), to_degrees(192, 0, 55), to_degrees(270, 12, 43),
               to_degrees(80, 57, 3), to_degrees(213, 6, 35)]
layings = [893.406, 514.887, 1146.720, 1664.159, 1460.275, 897.129, 836.043, 836.813]
first_point = [6065574.872, 4308709.395]
last_point = [6071450.856, 4307532.243]
first_dirangle = to_degrees(307, 7, 34)
last_dirangle = to_degrees(340, 1, 1)

