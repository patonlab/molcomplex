from os import path
import io
from setuptools import setup, find_packages

# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with io.open(path.join(this_directory, 'molcomplex', 'requirements.txt'), encoding='utf-8') as f:
    required = f.read().splitlines()

version = "1.0.1"
setup(
    name="molcomplex",
    packages=find_packages(),
    package_data={
        'molcomplex': ['metrics/*', 'metrics/*.gz'],
        'molcomplex': ['models/example_model/*','full_reaxys_model_1024bool/*',
                       'full_reaxys_model_1024uint8/*','full_reaxys_model_2048bool/*'],
    },
    version=version,
    license="MIT",
    description="Molecular Complexity Calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Guilian Luchini, Shree Sowndarya S. V.",
    author_email="guilian.luchini@colostate.edu, svss@colostate.edu",
    keywords=[
        "workflows",
        "computational chemistry",
        "cheminformatics",
        "molecular complexity",
        "DFT",
        "automation",
    ],
    url="https://github.com/patonlab/molcomplex",
    download_url=f"https://github.com/patonlab/molcomplex/archive/refs/tags/v1.0.1.tar.gz",
    install_requires=required,
    python_requires=">=3.0",
    include_package_data=True,
)
