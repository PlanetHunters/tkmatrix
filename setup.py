import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
version = "0.1.5"
setuptools.setup(
    name="tkmatrix", # Replace with your own username
    version=version,
    author="M. Dévora-Pajares & F.J. Pozuelos",
    author_email="mdevorapajares@protonmail.com",
    description="ToolKit for Multi-phase Augmented Transits Recovery from Injected eXoplanets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/martindevora/tkmatrix",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6.9',
    install_requires=['numpy==1.20.1; python_version>="3.7"',
                        'numpy==1.19; python_version<"3.7"',
                        "cython==0.29.21",
                        "pandas==1.1.5",
                        "lightkurve==2.0.2",
                        "wotan==1.9",
                        "matplotlib==3.3.4",
                        "corner==2.1.0",
                        "ellc==1.8.5",
                        "seaborn==0.11.1",
                        "sklearn==0.0",
                        "scipy==1.5.4",
                        "tess-point==0.6.1",
                        "astropy==4.1",
                        "mock==4.0.3",
                        'tqdm==4.56.0',
                        'setuptools>=41.0.0',
                        'torch==1.7.1',
                        'beautifulsoup4==4.9.3',
                        'numba>=0.53.0rc1',
                        'batman-package==2.4.7', # Transitleastsquares dependency
                        'argparse==1.4.0',
                        'configparser==5.0.1',
                        'pyparsing==2.4.7', # Matplotlib dependency
                        'transitleastsquares==1.0.25'
    ]
)