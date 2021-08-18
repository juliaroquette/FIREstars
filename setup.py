from setuptools import setup, find_packages


description = "Far-ultraviolet Irradiated Rotational Evolution"
description += "model for low mass stars"

setup(
    name='FIREstars',
    version='0.1',
    url="https://github.com/juliaroquette/FIREstars",
    author="Julia Roquette",
    author_email="juliaroquette@gmail.com",
    description=description,
    packages=find_packages(),
    include_package_data=True,
)
