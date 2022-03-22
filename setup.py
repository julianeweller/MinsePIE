import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="minsepie", 
    version="1.0",
    author="Juliane Weller",
    author_email="jw38@sanger.ac.uk",
    description="Tools to prime editing insertion efficiencies",
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
        'pandas',
        'scikit-learn>=0.24',
        'biopython>=1.79',
        'xgboost>=1.4',
        'datetime',
        'pandarallel>=1.6',
        'more-itertools>=8.12',
        'viennarna'
    ]
)