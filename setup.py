import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pycap",
    version="0.1.0",
    author="Michael N. Fienen, Aaron Pruitt, Howard Reeves",
    author_email="mnfienen@usgs.gov, aaron.pruitt@wisconsin.gov, hwreeves@usgs.gov",
    description="Stream Depletion analysis tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WDNR-Water-Use/HiCap_Analysis_Tool.",
    packages=setuptools.find_packages(),
    install_requires=["matplotlib","pyyaml","geopandas","numpy","scipy","pandas","gdal","fiona","rasterio>=1.0","rasterstats","shapely","rtree","pyproj>=2.0","pyshp","xlrd", "openpyxl", "requests", "pytest"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.11',
)
