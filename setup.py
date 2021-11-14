import setuptools

VERSION = "0.0.0"

setuptools.setup(
    name = "cytomulate",
    version = VERSION,
    description = "CyTOF Simulation",
    packages=["cytomulate"],
    python_requires=">=3.9",
    install_requires=["numpy",
                      "scipy"],
    test_requires=["pytest",
                   "pytest-cov",
                   "pytest-mock",
                   "coverage"],
    classifiers = [
        "Programming Language :: Python :: 3 :: Only",
        "Natural Language :: English"
    ]
)