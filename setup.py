from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

version = "1.1.0"
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
    long_description="Documentation in Read The Docs: XXX",
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
    download_url=f"https://github.com/patonlab/molcomplex/archive/refs/tags/v1.0.0.tar.gz",
    install_requires=required,
    python_requires=">=3.0",
    include_package_data=True,
)
