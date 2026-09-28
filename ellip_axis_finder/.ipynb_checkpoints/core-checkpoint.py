import numpy as np
from numpy.linalg import eigh
from scipy.stats import zscore
from sklearn.decomposition import PCA
from numba import njit
from joblib import Parallel, delayed

#----------------Method 1-----------------

@njit
def Compute_Distances_Numba(particle_pos, samp_center):
    """
    用Numba加速计算粒子位置与采样中心的欧氏距离
    
    参数:
        particle_pos: 形状为 (N, 3) 的NumPy数组，粒子的三维坐标
        samp_center: 形状为 (3,) 的NumPy数组，采样中心的三维坐标
    
    返回:
        distances: 形状为 (N,) 的NumPy数组，每个粒子到采样中心的距离
    """
    n = particle_pos.shape[0]
    distances = np.empty(n, dtype=np.float64)
    for i in range(n):
        dx = particle_pos[i, 0] - samp_center[0]
        dy = particle_pos[i, 1] - samp_center[1]
        dz = particle_pos[i, 2] - samp_center[2]
        distances[i] = np.sqrt(dx * dx + dy * dy + dz * dz)
    return distances


def Single_Sphere_Sampling(N_samp, sphere_center, radius, samp_radius, 
                          particle_pos, nearby_gal_pos, nearby_gal_rhalf, rho_N_ave, rho_N_cri_ratio):
    """
    单次球体采样并计算椭球轴（供并行调用）
    
    参数:
        N_samp: 目标采样点数
        sphere_center: 元组或形状为 (3,) 的数组，球体中心坐标 (x, y, z)
        radius: 球体半径
        samp_radius: 采样小球半径（= radius * samp_radius_ratio）
        particle_pos: 形状为 (N, 3) 的数组，粒子坐标
        rho_N_ave: 平均粒子数密度
        rho_N_cri_ratio: 密度阈值比例（实际阈值 = rho_N_cri_ratio * rho_N_ave）
    
    返回:
        axes: PCA主成分（3个轴方向），形状为 (3, 3)
        variances: 各主成分的方差（表示轴长度），形状为 (3,)
    """
    N_count = 0
    N_loop = 0
    samp_accept = []  # 存储接受的采样点 (x, y, z, 粒子数)
    
    while N_count < N_samp:
        # 解析球体中心坐标
        center_x, center_y, center_z = sphere_center
        
        # 生成球面上的随机采样点（球坐标系转笛卡尔坐标系）
        r = radius * np.random.random()  # 径向随机值（0~radius）
        theta = np.arccos(1 - 2 * np.random.random())  # 极角（0~π）
        phi = 2 * np.pi * np.random.random()  # 方位角（0~2π）
        
        # 转换为笛卡尔坐标
        x = center_x + r * np.sin(theta) * np.cos(phi)
        y = center_y + r * np.sin(theta) * np.sin(phi)
        z = center_z + r * np.cos(theta)
        samp_center = np.array([x, y, z])
        
        # 计算粒子到采样中心的距离并筛选
        distances = Compute_Distances_Numba(particle_pos, samp_center)
        samp_pos = particle_pos[distances < samp_radius]  # 筛选小球内的粒子
        
        # 密度判断：若密度超过阈值则接受该采样点
        volume = (4 / 3) * np.pi * (samp_radius **3)
        if volume > 0 and (len(samp_pos) / volume) >= rho_N_cri_ratio * rho_N_ave:
            samp_accept.append([samp_center[0], samp_center[1], samp_center[2], len(samp_pos)])
            N_count += 1
        
        N_loop = N_loop + 1
        if N_loop > 100 * N_samp:
            return None, None
    
    # 数据清洗（去除异常值）
    samp_accept = np.vstack(samp_accept)
    
    for i in range(len(nearby_gal_pos)):
        distances = Compute_Distances_Numba(samp_accept[:, 0:3], nearby_gal_pos[i, :])
        mask = distances >= 2 * nearby_gal_rhalf[i]
        samp_accept = samp_accept[mask]
        if len(samp_accept) < int(0.7 * N_samp):
            #print('The substructure is closely attached to the central galaxy. Skipping this halo!')
            return None, None
            
    z_scores = np.abs(zscore(samp_accept[:, 0:3]))  # 计算坐标的Z-score
    pos_clean = samp_accept[:, 0:3][(z_scores < 3).all(axis=1)]  # 保留Z-score < 3的点
    
    # PCA分析：计算椭球轴方向和方差
    pca = PCA(n_components=3)
    pca.fit(pos_clean)
    
    return pca.components_, pca.explained_variance_


def Sphere_Sampling_Method(N_samp, sphere_center, radius, samp_radius_ratio, 
                   particle_pos, nearby_gal_pos, nearby_gal_rhalf, rho_N_ave, rho_N_cri_ratio, N_re=1):
    """
    主函数：多进程并行球体采样，计算椭球轴分布
    
    参数:
        N_samp: 每次重复的目标采样点数
        sphere_center: 椭球中心坐标 (x, y, z)
        radius: 所有粒子分布空间的特征半径
        samp_radius_ratio: 采样小球半径比例（samp_radius = radius * 该值）
        particle_pos: (N, 3) 数组，粒子坐标
        rho_N_ave: 平均粒子数密度
        rho_N_cri_ratio: 密度阈值比例，采样小球内rho > rho_N_ave * rho_N_cri_ratio则接收
        N_re: 重复采样次数（并行执行）
    
    返回:
        axes_collect: 列表，每个元素为 (3, 3) 数组，存储每次重复的椭球三轴方向
        variances_collect: 列表，每个元素为 (3, N_re) 数组，存储每次重复的PCA方差
    """
    # 预处理参数
    samp_radius = radius * samp_radius_ratio
    particle_pos = np.ascontiguousarray(particle_pos)  # 确保内存连续，提升Numba效率
    
    # 并行执行N_re次采样
    results = Parallel(n_jobs=N_re, verbose=0)(  # 关闭verbose避免日志干扰
        delayed(Single_Sphere_Sampling)(
            N_samp=N_samp,
            sphere_center=sphere_center,
            radius=radius,
            samp_radius=samp_radius,
            particle_pos=particle_pos,
            nearby_gal_pos=nearby_gal_pos, 
            nearby_gal_rhalf=nearby_gal_rhalf, 
            rho_N_ave=rho_N_ave,
            rho_N_cri_ratio=rho_N_cri_ratio
        ) for _ in range(N_re)
    )
    
    for axes, variances in results:
        if axes is None or variances is None:
            return [None], [None]
        
    # 收集结果
    axes_collect = [res[0] for res in results]
    
    # 将第一个矢量作为标准，使重复计算的对应矢量与之夹角为锐角
    standard = axes_collect[0]
    adjusted_axes_collect = [standard.copy()]
    
    for i in range(1, len(axes_collect)):
        current_array = axes_collect[i].copy()
        # 检查每个矢量
        for j in range(current_array.shape[0]):
            # 计算当前矢量与标准矢量的点积
            dot_product = np.dot(current_array[j], standard[j])
            # 如果点积小于0，则反转当前矢量方向
            if dot_product < 0:
                current_array[j] *= -1
        adjusted_axes_collect.append(current_array)
                
    variances_collect = [res[1] for res in results]
    
    return adjusted_axes_collect, variances_collect


#----------------Method 2-----------------

@njit
def Compute_Inertia_Tensor(sphere_center, particle_pos, particle_mass, nearby_gal_pos, nearby_gal_rhalf): 
    for i in range(len(nearby_gal_pos)):
        distances = Compute_Distances_Numba(particle_pos[:, 0:3], nearby_gal_pos[i, :])
        mask = distances >= 2 * nearby_gal_rhalf[i]
        particle_pos = particle_pos[mask]
        particle_mass = particle_mass[mask]
    
    particle_pos = particle_pos - sphere_center
    inertia_tensor = np.zeros((3, 3), dtype=np.float64)
    
    # 遍历每个粒子
    for i in range(len(particle_mass)):
        mass = particle_mass[i]
        x, y, z = particle_pos[i]
        
        # 计算位置矢量的平方
        r_squared = x**2 + y**2 + z**2
        
        # 计算单个粒子对惯量张量的贡献并累加
        inertia_tensor[0, 0] += mass * (r_squared - x**2)
        inertia_tensor[0, 1] += mass * (-x * y)
        inertia_tensor[0, 2] += mass * (-x * z)
        
        inertia_tensor[1, 0] += mass * (-y * x)
        inertia_tensor[1, 1] += mass * (r_squared - y**2)
        inertia_tensor[1, 2] += mass * (-y * z)
        
        inertia_tensor[2, 0] += mass * (-z * x)
        inertia_tensor[2, 1] += mass * (-z * y)
        inertia_tensor[2, 2] += mass * (r_squared - z**2)
    
    return inertia_tensor

def Inertia_Tensor_Method(sphere_center, particle_pos, particle_mass, nearby_gal_pos, nearby_gal_rhalf):
    inertia_tensor = Compute_Inertia_Tensor(sphere_center, particle_pos, particle_mass, nearby_gal_pos, nearby_gal_rhalf)
    if inertia_tensor is None:
        return [None], [None]
    moments, axes = eigh(inertia_tensor)
    return axes.T, moments