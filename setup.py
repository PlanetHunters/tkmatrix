import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
version = "0.15.1"

setuptools.setup(
    name="tkmatrix",
    version=version,
    author="M. Dévora-Pajares & F.J. Pozuelos",
    author_email="mdevorapajares@protonmail.com",
    description="ToolKit for Multi-phAse Transits Recovery from Injected eXoplanets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PlanetHunters/tkmatrix",
    license_files=["LICENSE", "tkmatrix/ravex/LICENSE"],
    # ravex is a git submodule keeping its code under src/, which find_packages does not reach
    packages=setuptools.find_packages() + ['tkmatrix.ravex.src.ravex'],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.11',
    install_requires=['argparse==1.4.0',
                        'beautifulsoup4==4.13.4',
                        'configparser==5.0.1',
                        "corner==2.2.3",
                        "ellc==1.9.0", # Light curve models (Fortran, built with gfortran at install time)
                        "lcbuilder==0.26.0",
                        "mock==5.2.0",
                        'pyparsing==3.2.3', # Matplotlib dependency
                        "seaborn==0.13.2",
                        'setuptools>=41.0.0',
                        "sklearn==0.0"
    ]
)
