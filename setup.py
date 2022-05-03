from setuptools import setup

VERSION = "0.1.0"

setup(
    name = "cytomulate",
    version = VERSION,
    description = "A package for simulating CyTOF data",
    author = "Yuqiu Yang, Kevin Wang, Tao Wang, Sherry Wang",
    author_email = "yuqiuy@smu.edu, kevinwang@smu.edu, Tao.Wang@UTSouthwestern.edu, swang@mail.smu.edu",
    packages=["cytomulate"],
    python_requires=">=3.5",
    install_requires=["numpy",
                      "scipy",
                      "scikit-learn",
                      "networkx",
                      "matplotlib",
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
)