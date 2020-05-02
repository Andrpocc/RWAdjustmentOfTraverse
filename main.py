from parametric_adjust import parametric_adjustment, adjustment_traverse, to_degrees
import numpy as np
import pandas as pd
from random import normalvariate
from math import sqrt

input_angles = [to_degrees(103, 44, 10), to_degrees(209, 43, 22), to_degrees(180, 14, 38),
                to_degrees(127, 38, 47), to_degrees(182, 13, 8), to_degrees(202, 30, 15), to_degrees(182, 50, 32),
                to_degrees(216, 47, 44), to_degrees(96, 10, 43)]
input_layings = [410.432, 435.832, 400.294, 474.048, 500.093, 373.122, 420.237, 605.386]
input_first_point = [6066110.881, 4310283.192]
input_last_point = [6069429.370, 4310505.674]
first_dirangle = to_degrees(71, 5, 8)
last_dirangle = to_degrees(312, 58, 27)
input_x_coordinates = [6066519.638, 6066916.088, 6067279.500, 6067699.735, 6068151.687, 6068524.350, 6068944.587]
input_y_coordinates = [4310246.148, 4310427.193, 4310595.023, 4310375.654, 4310161.569, 4310143.068, 4310143.068]
input_coordinates = pd.DataFrame()
input_coordinates['x'] = input_x_coordinates
input_coordinates['y'] = input_y_coordinates

x_matrix = pd.DataFrame()
y_matrix = pd.DataFrame()

for i in range(2):
    mistakes = list()
    for j in range(9):
        mistakes.append(normalvariate(0, 5) / 3600)
    mistake_angles = np.array(input_angles) + np.array(mistakes)
    mistakes = list()
    for j in range(8):
        mistakes.append(normalvariate(0, 2 + 2 * input_layings[j] / 1000) / 1000)
    mistake_layings = np.array(input_layings) + np.array(mistakes)
    corrected_coordinates = adjustment_traverse(first_dirangle, last_dirangle, mistake_angles, mistake_layings,
                                                input_first_point, input_last_point)
    delta_x = np.array(input_coordinates['x'] - corrected_coordinates['x'])
    delta_y = np.array(input_coordinates['y'] - corrected_coordinates['y'])
    delta_x_kv = np.power(delta_x, 2)
    delta_y_kv = np.power(delta_y, 2)
    x_matrix[str(i)] = delta_x_kv
    y_matrix[str(i)] = delta_y_kv

mt_list = list()
mx_list = list()
my_list = list()
for i in range(7):
    mx = sqrt(sum(x_matrix.values[i]) / len(x_matrix.values[i]))
    mx_list.append(mx)
for i in range(7):
    my = sqrt(sum(y_matrix.values[i]) / len(y_matrix.values[i]))
    my_list.append(my)
for i in range(7):
    mt = sqrt(mx_list[i] ** 2 + my_list[i] ** 2)
    mt_list.append(mt)

print(mt_list)
