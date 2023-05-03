from setuptools import setup,find_packages
setup(
    name="NanoSPA",
    version="1.0",
    description="Nanopore m6A and pseudouridine sequencing",
    author="Sihao Huang",
    packages=find_packages(),
    include_package_data=True,
    scripts=["bin/nanospa"],
    license="GPL 2.0"
)