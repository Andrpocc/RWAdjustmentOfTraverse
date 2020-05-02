import pandas as pd
import numpy as np
from random import normalvariate
from geo_tasks import *
from math import sqrt, sin, cos, radians, degrees, atan2


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
    theoretical_angles.append(ogz_points(first_point[0], first_point[1], x_coordinates[0],
                                         y_coordinates[0]) - first_directional_angle + 180 - 360)
    theoretical_angles.append(
        degrees(atan2(y_coordinates[1] - y_coordinates[0], x_coordinates[1] - x_coordinates[0])) - degrees(
            atan2(first_point[1] - y_coordinates[0], first_point[0] - x_coordinates[0])) + 360)
    for i in range(len(angles_array) - 4):
        angle = degrees(
            atan2(y_coordinates[i + 2] - y_coordinates[i + 1], x_coordinates[i + 2] - x_coordinates[i + 1])) - degrees(
            atan2(y_coordinates[i] - y_coordinates[i + 1], x_coordinates[i] - x_coordinates[i + 1]))
        if angle < 0:
            theoretical_angles.append(angle + 360)
        elif angle > 360:
            theoretical_angles.append(angle - 360)
        else:
            theoretical_angles.append(angle)

    theoretical_angles.append(
        degrees(atan2(last_point[1] - y_coordinates[-1], last_point[0] - x_coordinates[-1])) - degrees(
            atan2(y_coordinates[-2] - y_coordinates[-1], x_coordinates[-2] - x_coordinates[-1])) + 360)
    theoretical_angles.append(last_directional_angle - ogz_points(x_coordinates[-1], y_coordinates[-1], last_point[0],
                                                                  last_point[1]) + 180 - 360)
    print(theoretical_angles)
    print(angles)
    # Теоретические значения горизонтальных проложений
    theoretical_layings = list()
    theoretical_layings.append(
        sqrt((x_coordinates[0] - first_point[0]) ** 2 + (y_coordinates[0] - first_point[1]) ** 2))
    for i in range(len(horizontal_layings_array) - 2):
        theoretical_layings.append(
            sqrt((x_coordinates[i + 1] - x_coordinates[i]) ** 2 + (y_coordinates[i + 1] - y_coordinates[i]) ** 2))
    theoretical_layings.append(
        sqrt((last_point[0] - x_coordinates[-1]) ** 2 + (last_point[1] - y_coordinates[-1]) ** 2))
    # Составление матрицы свободных членов L
    v_angles = (np.array(theoretical_angles) - np.array(angles)) * 3600
    v_layings = (np.array(theoretical_layings) - np.array(horizontal_layings)) * 100
    l_matrix = np.concatenate((v_angles, v_layings))  # Поменял знак, так как перенос на другую сторону
    # Составлание матрицы коэффициентов А для углов
    a_matrix_angles = list()
    for i in range(len(v_angles)):
        row = list()
        for j in range(len(x_coordinates) * 2):
            row.append(0)
        a_matrix_angles.append(row)
    # Коэффициенты для поправок в углы
    # Коэффиценты для первого угла
    kxj = -(sin(radians(directional_angles[0])) * (206265 / (horizontal_layings[0] * 100)))
    kyj = (cos(radians(directional_angles[0])) * (206265 / (horizontal_layings[0] * 100)))
    a_matrix_angles[0][0] = kxj
    a_matrix_angles[0][1] = kyj
    # Теперь цикл для всех углов до последнего
    for i in range(len(v_angles) - 2):
        if directional_angles[i] < 180:
            dir_angle_ik = directional_angles[i] + 180
        else:
            dir_angle_ik = directional_angles[i] - 180
        kxk = (sin(radians(dir_angle_ik)) * (206265 / (horizontal_layings[i] * 100)))
        kyk = -(cos(radians(dir_angle_ik)) * (206265 / (horizontal_layings[i] * 100)))
        kxi = (sin(radians(directional_angles[i + 1])) * (206265 / (horizontal_layings[i+1] * 100))) - (
                sin(radians(dir_angle_ik)) * (206265 / (horizontal_layings[0] * 100)))
        kyi = -(cos(radians(directional_angles[i + 1])) * (206265 / (horizontal_layings[i+1] * 100))) + (
                cos(radians(dir_angle_ik)) * (206265 / (horizontal_layings[0] * 100)))
        kxj = -(sin(radians(directional_angles[i + 1])) * (206265 / (horizontal_layings[i+1] * 100)))
        kyj = (cos(radians(directional_angles[i + 1])) * (206265 / (horizontal_layings[i+1] * 100)))
        # Заполение матрицы
        if i == 0:
            a_matrix_angles[1][0] = kxi
            a_matrix_angles[1][1] = kyi
            a_matrix_angles[1][2] = kxj
            a_matrix_angles[1][3] = kyj
            point = 0  # Для перехода точек через две координаты
        elif i == (len(v_angles) - 3):
            a_matrix_angles[i + 1][point] = kxk
            a_matrix_angles[i + 1][point + 1] = kyk
            a_matrix_angles[i + 1][point + 2] = kxi
            a_matrix_angles[i + 1][point + 3] = kyi
        else:
            a_matrix_angles[i + 1][point] = kxk
            a_matrix_angles[i + 1][point + 1] = kyk
            a_matrix_angles[i + 1][point + 2] = kxi
            a_matrix_angles[i + 1][point + 3] = kyi
            a_matrix_angles[i + 1][point + 4] = kxj
            a_matrix_angles[i + 1][point + 5] = kyj
            point += 2
    # Коэффиценты последнего угла
    if directional_angles[-1] < 180:
        dir_angle_ik = directional_angles[-1] + 180
    else:
        dir_angle_ik = directional_angles[-1] - 180
    kxk = (sin(radians(dir_angle_ik)) * (206265 / (horizontal_layings[-1] * 100)))
    kyk = -(cos(radians(dir_angle_ik)) * (206265 / (horizontal_layings[-1] * 100)))
    a_matrix_angles[-1][-2] = kxk
    a_matrix_angles[-1][-1] = kyk
    # Составление матрицы А для проложений
    a_matrix_layings = list()
    for i in range(len(v_layings)):
        row = list()
        for j in range(len(x_coordinates) * 2):
            row.append(0)
        a_matrix_layings.append(row)
    # Коэфициенты для поправок в стороны
    # Коэфициенты для первой стороны
    kxj = (cos(radians(directional_angles[0])))
    kyj = (sin(radians(directional_angles[0])))
    a_matrix_layings[0][0] = kxj
    a_matrix_layings[0][1] = kyj
    # Теперь цикл для всех сторон до последнего
    for i in range(len(v_layings) - 2):
        kxi = -(cos(radians(directional_angles[i + 1])))
        kyi = -(sin(radians(directional_angles[i + 1])))
        kxj = cos(radians(directional_angles[i + 1]))
        kyj = sin(radians(directional_angles[i + 1]))
        # Заполение матрицы
        if i == 0:
            point = 0  # Для перехода точек через две координаты
        a_matrix_layings[i + 1][point] = kxi
        a_matrix_layings[i + 1][point + 1] = kyi
        a_matrix_layings[i + 1][point + 2] = kxj
        a_matrix_layings[i + 1][point + 3] = kyj
        point += 2
    # Коэффициенты последей стороны
    kxi = -(cos(radians(directional_angles[-1])))
    kyi = -(sin(radians(directional_angles[-1])))
    a_matrix_layings[-1][-2] = kxi
    a_matrix_layings[-1][-1] = kyi
    # Построение матрицы общей матрицы А из матрицы коэф углои и матрицы коэф сторон
    a_matrix_angles_array = np.array(a_matrix_angles)
    a_matrix_layings_array = np.array(a_matrix_layings)
    a_matrix_array = np.concatenate((a_matrix_angles_array, a_matrix_layings_array))
    # Составление матрицы весов
    p_matrix = list()
    for i in range(len(l_matrix)):
        row = list()
        for j in range(len(l_matrix)):
            row.append(0)
        p_matrix.append(row)
    # Угловые веса
    for i in range(len(v_angles)):
        p_matrix[i][i] = 1
    # Линейные веса
    for i in range(len(v_angles), len(l_matrix)):
        p_matrix[i][i] = 5 ** 2 / 0.2 ** 2
    p_matrix_array = np.array(p_matrix)
    # Составление матрицы нормальных уравнений RA
    a_matrix_tr_array = np.transpose(a_matrix_array)
    r_koeff = np.dot(a_matrix_tr_array, p_matrix_array)
    ra = np.dot(r_koeff, a_matrix_array)
    rl = -np.dot(r_koeff, l_matrix)
    # Решение матричноего уравнения поправки в координаты
    v_coordinates = np.linalg.solve(ra, rl)
    # Нахождение поправок в измерения
    v_measurement = np.dot(a_matrix_array, v_coordinates) + l_matrix
    # контроль
    control = np.dot(np.dot(a_matrix_tr_array, p_matrix_array), v_measurement)
    pvv = np.dot(np.dot(np.transpose(v_measurement), p_matrix_array), v_measurement)
    pvl = np.dot(np.dot(np.transpose(v_measurement), p_matrix_array), l_matrix)
    # Вычисление уравненных коодинат
    v_coordinates = list(v_coordinates)
    v_x = list()
    for i in range(0, len(v_coordinates), 2):
        v_x.append(v_coordinates[i])
    v_y = list()
    for i in range(1, len(v_coordinates), 2):
        v_y.append(v_coordinates[i])
    v_x = np.array(v_x)
    v_y = np.array(v_y)
    x_coordinates = np.array(x_coordinates)
    y_coordinates = np.array(y_coordinates)
    x_corrected = x_coordinates + v_x / 100
    y_corrected = y_coordinates + v_y / 100
    table = pd.DataFrame()
    table['x'] = x_corrected
    table['y'] = y_corrected
    # Оценка точности
    mu = sqrt(pvv / 3)
    mb = mu
    ms = mu * sqrt(1 / (10 ** 2 / 2 ** 2))
    k = np.linalg.inv(ra) * mu ** 2
    q = list()
    for i in range(len(v_coordinates)):
        q.append(sqrt(k[i][i]))
    mt = list()
    for i in range(0, len(q), 2):
        mt.append(sqrt(q[i] ** 2 + q[i + 1] ** 2))
    print(v_coordinates)
    print(l_matrix)
    print(mt)
    return l_matrix


if __name__ == '__main__':
    left_angles = [to_degrees(103, 44, 13), to_degrees(209, 43, 14), to_degrees(180, 14, 35),
                   to_degrees(127, 38, 45), to_degrees(182, 13, 15), to_degrees(202, 29, 58), to_degrees(182, 50, 33),
                   to_degrees(216, 47, 46), to_degrees(96, 10, 44)]
    layings = [410.438, 435.831, 400.290, 474.049, 500.093, 373.122, 420.239, 605.389]
    first_point = [6066110.881, 4310283.192]
    last_point = [6069429.370, 4310505.674]
    first_dirangle = to_degrees(71, 5, 8)
    last_dirangle = to_degrees(312, 58, 27)

    mistakes1 = list()
    for j in range(9):
        mistakes1.append(normalvariate(0, 5))
    mistakes2 = list()
    for j in range(8):
        mistakes2.append(normalvariate(0, 2 + 2 * layings[j] / 1000))
    with open('mistakes.txt', 'w') as file:
        print(mistakes1, file=file)
        print(mistakes2, file=file)
    table = parametric_adjustment(first_dirangle, last_dirangle, left_angles, layings, first_point,
                                  last_point, left_angle=True)
