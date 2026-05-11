import setuptools


setuptools.setup(
    name="YAOBridge",
    version="1.1.0",
    author="kangLab yijiyang",
    author_email="15623109189@163.com",
    description="A chain tool for confident position lift between hg38 and yao",
    url="",
    install_requires=[],
    python_requires='>=3.10',
    packages=setuptools.find_packages(),
include_package_data=True, # 建议开启
    package_data={
        # 表示在 src 包下的 data 目录中包含所有 .chain 文件
        "src": ["data/*.chain"], 
    },
    entry_points={'console_scripts': ['YAOBridge = src.PB:main'], },
)
