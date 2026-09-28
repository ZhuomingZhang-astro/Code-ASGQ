import numpy as np
import h5py

def SFGVQ(sample_angle, sgal_arr, SFgal, GVgal, Qgal, snap_num, redshift, bins_num, snapnum_center, base_save_path, save=True):
    # 【进一步】筛选 Starforming, Green Valley 与 Quenched星系：   
    for i in range(len(sgal_arr)):
        if (sgal_arr[i][6] != 0) and (np.log10(sgal_arr[i][6]) > (0.73 * np.log10(sgal_arr[i][7])\
                                                                                  - 7.33 + np.log10((1 + redshift)**2))):
            SFgal.append(sgal_arr[i][:])

        if (sgal_arr[i][6] != 0) and (np.log10(sgal_arr[i][6]) < (0.73 * np.log10(sgal_arr[i][7])\
                                                                                  - 7.33 + np.log10((1 + redshift)**2))) and \
        (np.log10(sgal_arr[i][6]) > (0.73 * np.log10(sgal_arr[i][7]) - 7.33 - 1 + np.log10((1 + redshift)**2))):

            GVgal.append(sgal_arr[i][:])    

        if (sgal_arr[i][6] == 0) or (np.log10(sgal_arr[i][6]) < (0.73 * np.log10(sgal_arr[i][7])\
                                                                                 - 7.33 - 1 + np.log10((1 + redshift)**2))):
            Qgal.append(sgal_arr[i][:])

    if SFgal:
        SFgal_arr = np.vstack(SFgal)
    else:
        SFgal_arr = np.array([])  # or handle it as needed
    
    
    if GVgal:
        GVgal_arr = np.vstack(GVgal)
    else:
        GVgal_arr = np.array([])  # or handle it as needed

    if Qgal:
        Qgal_arr = np.vstack(Qgal)
    else:
        Qgal_arr = np.array([])


    orient_angle = 30
    orient = np.linspace(0, 360, int((360/orient_angle)+1))
    x = np.cos(np.radians(orient))
    y = np.zeros(len(orient))
    z = np.sin(np.radians(orient))
    orient_vec= np.column_stack((x, y, z))

    # 初始化一个列表来存储每个圆锥内的粒子
    SFinside_cone = [[] for _ in range(len(orient_vec))]
    GVinside_cone = [[] for _ in range(len(orient_vec))]
    Qinside_cone = [[] for _ in range(len(orient_vec))]

    # 遍历每个圆锥
    for i in range(len(orient_vec)):
        for j in range(SFgal_arr.shape[0]):
            # 计算粒子与圆锥轴向量的点积
            dot_products = np.dot(SFgal_arr[j, 0:3], orient_vec[i])
            # 计算角度（以弧度为单位）
            angles = np.arccos(dot_products / (np.linalg.norm(orient_vec[i]) * np.linalg.norm(SFgal_arr[j, 0:3,])))
            # 将角度转换为度
            angles_deg = np.degrees(angles)
            if (angles_deg < sample_angle):
                SFinside_cone[i].append(SFgal_arr[j, :])
        # 检查是否为空
        if SFinside_cone[i]:
            SFinside_cone[i] = np.vstack(SFinside_cone[i])
        else:
            SFinside_cone[i] = np.empty((0, sgal_arr[0].shape[0]))  # 或者使用其他适当的值

    for i in range(len(orient_vec)):
        for j in range(GVgal_arr.shape[0]):
            # 计算粒子与圆锥轴向量的点积
            dot_products = np.dot(GVgal_arr[j, 0:3], orient_vec[i])
            # 计算角度（以弧度为单位）
            angles = np.arccos(dot_products / (np.linalg.norm(orient_vec[i]) * np.linalg.norm(GVgal_arr[j, 0:3])))
            # 将角度转换为度
            angles_deg = np.degrees(angles)
            if (angles_deg < sample_angle):
                GVinside_cone[i].append(GVgal_arr[j, :])
        # 检查是否为空
        if GVinside_cone[i]:
            GVinside_cone[i] = np.vstack(GVinside_cone[i])
        else:
            GVinside_cone[i] = np.empty((0, sgal_arr[0].shape[0]))  # 或者使用其他适当的值

    for i in range(len(orient_vec)):
        for j in range(Qgal_arr.shape[0]):
            # 计算粒子与圆锥轴向量的点积
            dot_products = np.dot(Qgal_arr[j, 0:3], orient_vec[i])
            # 计算角度（以弧度为单位）
            angles = np.arccos(dot_products / (np.linalg.norm(orient_vec[i]) * np.linalg.norm(Qgal_arr[j, 0:3])))
            # 将角度转换为度
            angles_deg = np.degrees(angles)
            if (angles_deg < sample_angle):
                Qinside_cone[i].append(Qgal_arr[j, :])
        # 检查是否为空
        if Qinside_cone[i]:
            Qinside_cone[i] = np.vstack(Qinside_cone[i])
        else:
            Qinside_cone[i] = np.empty((0, sgal_arr[0].shape[0]))  # 或者使用其他适当的值


    cone_data = np.empty(3, dtype=object)
    cone_data[0] = SFinside_cone
    cone_data[1] = GVinside_cone
    cone_data[2] = Qinside_cone
    
    np.save(base_save_path + f'/cone_{bins_num}_{snap_num:03}.npy', cone_data)



    Ndata = []
    for i in range(len(orient_vec)):
        Ndata.append((SFinside_cone[i].shape[0] + GVinside_cone[i].shape[0] + Qinside_cone[i].shape[0]))
    Ndata = np.array(Ndata)

    QNdata = []
    for i in range(len(orient_vec)):
        QNdata.append((Qinside_cone[i].shape[0]))
    QNdata = np.array(QNdata)

    GVNdata = []
    for i in range(len(orient_vec)):
        GVNdata.append((GVinside_cone[i].shape[0]))
    GVNdata = np.array(GVNdata)

    SFNdata = []
    for i in range(len(orient_vec)):
        SFNdata.append((SFinside_cone[i].shape[0]))
    SFNdata = np.array(SFNdata)

    Qfracdata = []
    for i in range(len(orient_vec)):
        if ((SFinside_cone[i].shape[0] + GVinside_cone[i].shape[0] + Qinside_cone[i].shape[0]) != 0):
            Qfracdata.append(Qinside_cone[i].shape[0] / (SFinside_cone[i].shape[0] + GVinside_cone[i].shape[0] + Qinside_cone[i].shape[0]))
        else:
            Qfracdata.append(np.nan)
    Qfracdata = np.array(Qfracdata)


    N_distri = np.vstack((Ndata, QNdata, GVNdata, SFNdata, Qfracdata, orient))
    np.save(base_save_path + f'/NQfrac_{bins_num}_{snap_num:03}.npy', N_distri)
        
    return