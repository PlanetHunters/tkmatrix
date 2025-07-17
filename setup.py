import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
version = "0.12.0"
import subprocess
import shutil
import os
from setuptools.command.install import install

class CustomInstall(install):
    def run(self):
        here = os.path.abspath(os.path.dirname(__file__))
        ellc_dir = os.path.join(here, 'tkmatrix', 'ellc')
        print(f"=== Running make in {ellc_dir} ===")
        subprocess.check_call(['make', '-B'], cwd=ellc_dir)
        shutil.copy(ellc_dir + '/libellc.so', os.path.join(ellc_dir, 'ellc') + '/libellc.so')
        super().run()

setuptools.setup(
    name="tkmatrix",
    version=version,
    author="M. Dévora-Pajares & F.J. Pozuelos",
    author_email="mdevorapajares@protonmail.com",
    description="ToolKit for Multi-phAse Transits Recovery from Injected eXoplanets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PlanetHunters/tkmatrix",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    cmdclass={'install': CustomInstall},
    python_requires='>=3.11',
    install_requires=['argparse==1.4.0',
                        'beautifulsoup4==4.13.4',
                        'configparser==5.0.1',
                        "corner==2.2.3",
                        "lcbuilder==0.25.4",
                        "mock==5.2.0",
                        'pyparsing==3.2.3', # Matplotlib dependency
                        "seaborn==0.13.2",
                        'setuptools>=41.0.0',
                        "sklearn==0.0"
    ]
)
