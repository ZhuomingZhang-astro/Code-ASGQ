from setuptools import setup, find_packages

setup(
    name="correct_periodic_coords",  # 包名
    version="0.1.1",  # 版本号
    packages=find_packages(),  # 自动发现包内所有模块
    install_requires=[  # 依赖的第三方库
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "numba>=0.55.0",
        "joblib>=1.1.0"
    ],
    author="Zhuoming Zhang",  # 作者名（可修改）
    description="A tool for finding ellipsoid axes via sphere sampling and PCA",  # 包描述
    url="",  # 项目地址（可选）
)