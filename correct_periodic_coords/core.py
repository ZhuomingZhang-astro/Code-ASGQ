import numpy as np
from numba import njit

@njit
def Correct_Periodic_Coords(pos, central_pos, radius, L_box):

    # 计算参数
    h = 0.68
    threshold = 30 * radius
    compensate_amp = L_box
    
    # 创建结果数组（避免修改原数组）
    result = np.empty_like(pos)
    
    # 显式循环处理（Numba对循环优化极好）
    for i in range(len(pos)):
        for d in range(len(central_pos)):
            delta = pos[i, d] - central_pos[d]
            if np.abs(delta) > threshold:
                # 计算符号并应用补偿
                sign = 1.0 if delta > 0 else -1.0
                result[i, d] = pos[i, d] - sign * compensate_amp
            else:
                result[i, d] = pos[i, d]
    
    # 应用补偿，返回修正后的坐标
    return result