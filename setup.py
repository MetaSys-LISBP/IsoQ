import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="IsoQ",
    version="0.1",
    author="Pierre Millard, Baudoin Delépine",
    author_email="millard@insa-toulouse.fr",
    description="Description.....",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MetaSys-LISBP/IsoQ/",
    packages=setuptools.find_packages(),
    python_requires='>=3.5',
    install_requires=['pandas>=0.17.1'],
    package_data={'': ['data/*.csv', ], },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ]
)
