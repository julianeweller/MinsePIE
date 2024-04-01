import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="minsepie", 
    version="3.1",
    author="Juliane Weller",
    author_email="jw38@sanger.ac.uk",
    description="Tool to prime editing insertion efficiencies",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'pandas>=1.3',
        'regex>=2021.8',
        'scikit-learn>=0.24',
        'biopython>=1.79',
        'xgboost==1.5.0',
        'scipy>=1.5',
        'datetime',
        'pandarallel==1.5.4',
        'more-itertools>=8.12',
        'viennarna>=2.5.0a1',
        'psutil==5.9.0',
        'more_itertools>=8.12.0'
    ]
)
