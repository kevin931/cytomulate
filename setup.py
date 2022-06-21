from setuptools import setup
import sys
import os
import shutil
import distutils.cmd

VERSION = "0.0.1"

class PypiCommand(distutils.cmd.Command):
    
    description = "Build and upload for PyPI."
    user_options = []
    
    def initialize_options(self):
        pass
    
    
    def finalize_options(self):
        pass
    
    
    def run(self):
        try:
            shutil.rmtree("dist/")
        except FileNotFoundError:
            pass
        
        wheel_file = "cytomulate-{}-py3-none-any.whl".format(VERSION)
        tar_file = "cytomulate-{}.tar.gz".format(VERSION)
        
        os.system("{} setup.py sdist bdist_wheel".format(sys.executable))
        os.system("twine upload dist/{} dist/{}".format(wheel_file, tar_file))
    
    
class CondaCommand(distutils.cmd.Command):
    
    description = "Build and upload for conda."
    user_options = []
    
    def initialize_options(self):
        pass
    
    
    def finalize_options(self):
        pass
    
    
    def run(self):
        try:
            shutil.rmtree("dist_conda/")
        except FileNotFoundError:
            pass
        os.system("conda build . --output-folder dist_conda/")
        os.system("anaconda upload ./dist_conda/noarch/cytomulate-{}-py_0.tar.bz2".format(VERSION))


setup(
    name = "cytomulate",
    version = VERSION,
    description = "Accurate and Efficient Simulation of CyTOF data",
    author = "Yuqiu Yang, Kevin Wang, Tao Wang, Sherry Wang",
    author_email = "yuqiuy@smu.edu, kevinwang@smu.edu, Tao.Wang@UTSouthwestern.edu, swang@mail.smu.edu",
    long_description_content_type = "text/markdown",
    long_description = open("README.md").read(),
    packages=["cytomulate"],
    python_requires=">=3.5",
    install_requires=["numpy",
                      "scipy",
                      "scikit-learn",
                      "networkx",
                      "matplotlib",
                      "tqdm"
                      ],
    test_requires=["pytest",
                   "pytest-cov",
                   "pytest-mock",
                   "coverage",
                   ],
    classifiers = [
        "Programming Language :: Python :: 3 :: Only",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        ],
    cmdclass = {"pypi": PypiCommand,
                "conda": CondaCommand
                }
)