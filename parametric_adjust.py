import pandas as pd
import numpy as np
from geo_tasks import *
import pprint
from math import sqrt

def parametric_adjustment(first_directional_angle: float, last_directional_angle: float, angles: list,
                          horizontal_layings: list, first_point: list, last_point: list, left_angle=True):
    # Вычисление значений приближенных координат пунктов полигонометрического хода
    angles_array = np.array(angles)
    horizontal_layings_array = np.array(horizontal_layings)
    if left_angle:
        directional_angles = list()
        directional_angles.append(first_directional_angle + angles_array[0] - 180)
        for i in range(len(angles_array) - 2):
            dir_angle = directional_angles[i] + angles_array[i + 1] - 180
            if dir_angle > 360:
                directional_angles.append(dir_angle - 360)
            elif dir_angle < 0:
                directional_angles.append(dir_angle + 360)
            else:
                directional_angles.append(dir_angle)
    else:
        directional_angles = list()
        directional_angles.append(first_directional_angle - angles_array[0] + 180)
        for i in range(len(angles_array) - 2):
            directional_angles.append(directional_angles[i] - angles_array[i + 1] + 180)
    directional_angles_array = np.array(directional_angles)
    delta_x_array = horizontal_layings_array * np.cos(np.radians(directional_angles_array))
    delta_y_array = horizontal_layings_array * np.sin(np.radians(directional_angles_array))
    x_coordinates = list()
    x_coordinates.append(round(first_point[0] + delta_x_array[0], 3))
    for i in range(len(delta_x_array) - 2):
        x_coordinates.append(x_coordinates[i] + delta_x_array[i + 1])
    y_coordinates = list()
    y_coordinates.append(round(first_point[1] + delta_y_array[0], 3))
    for i in range(len(delta_y_array) - 2):
        y_coordinates.append(y_coordinates[i] + delta_y_array[i + 1])
    # По полученным приближенным значениям кооодинат вычисляются теоретические углы и длины линий
    # Теоретические значения горизонтальных углов
    theoretical_angles = list()
    theoretical_angles.append(directional_angles[0] - first_directional_angle + 180)
    for i in range(len(angles_array) - 2):
        angle = directional_angles[i + 1] - directional_angles[i] + 180
        if angle < 0:
            theoretical_angles.append(angle + 360)
        elif angle > 360:
            theoretical_angles.append(angle - 360)
        else:
            theoretical_angles.append(angle)
    theoretical_angles.append(last_directional_angle - directional_angles[-1] + 180)
    # Теоретические значения горизонтальных проложений
    theoretical_layings = list()
    theoretical_layings.append(sqrt((x_coordinates[0] - first_point[0])**2 + (y_coordinates[0] - first_point[1])**2))
    for i in range(len(horizontal_layings_array) - 2):
        theoretical_layings.append(
            sqrt((x_coordinates[i + 1] - x_coordinates[i]) ** 2 + (y_coordinates[i + 1] - y_coordinates[i]) ** 2))
    theoretical_layings.append(sqrt((last_point[0] - x_coordinates[-1])**2 + (last_point[1] - y_coordinates[-1])**2))
    # Составление матрицы свободных членов L
    v_angles = np.array(theoretical_angles) - np.array(angles)
    v_layings = np.array(theoretical_layings) - np.array(horizontal_layings)
    l_matrix = -np.concatenate((v_angles, v_layings))  # Поменял знак, так как перенос на другую сторону
    print(v_angles)
    print(l_matrix)
    # Составлание матрицы коэффициентов А
    # Коэффициенты для поправок в углы



    table =pd.DataFrame()
    table['angl'] = theoretical_angles
    table['la'] = angles
    #print(table)



if __name__ == '__main__':
    left_angles = [to_degrees(206, 15, 10), to_degrees(146, 2, 40), to_degrees(285, 26, 55),
                   to_degrees(162, 44, 24), to_degrees(96, 6, 53), to_degrees(192, 0, 51), to_degrees(270, 12, 41),
                   to_degrees(80, 57, 15), to_degrees(213, 6, 39)]
    layings = [893.4, 514.874, 1146.725, 1664.18, 1460.28, 897.177, 835.992, 836.817]
    first_point = [6065574.872, 4308709.395]
    last_point = [6071450.856, 4307532.243]
    parametric_adjustment(to_degrees(307, 7, 34), to_degrees(340, 1, 1), left_angles, layings, first_point,
                          last_point, left_angle=True)
