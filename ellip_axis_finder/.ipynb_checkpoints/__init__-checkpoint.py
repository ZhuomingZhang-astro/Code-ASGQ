"""
ellip_axis_finder: 球体采样与椭球轴计算工具包
可用通过[球体随机采样和PCA分析]或[惯量主轴计算]计算粒子分布的椭球轴方向和长度
"""

from .core import (
    Compute_Distances_Numba,
    Single_Sphere_Sampling,
    Sphere_Sampling_Method,
    Compute_Inertia_Tensor,
    Inertia_Tensor_Method
)

__all__ = [
    "Compute_Distances_Numba",
    "Single_Sphere_Sampling",
    "Sphere_Sampling_Method",
    "Compute_Inertia_Tensor",
    "Inertia_Tensor_Method"
]

__version__ = "0.5.7"  # 包版本号