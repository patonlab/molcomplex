from os import path
import io
from setuptools import setup, find_packages

# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

version = "1.0.3"
setup(
    name="molcomplex",
    packages=find_packages(),
    package_data={
        'molcomplex': ['requirements.txt'],
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
    download_url=f"https://github.com/patonlab/molcomplex/archive/refs/tags/v1.0.3.tar.gz",
    install_requires=[
        'brotli==1.1.0',
        'click==8.1.3',
        'dash-bootstrap-components==1.6.0',
        'dash-core-components==2.0.0',
        'dash-dangerously-set-inner-html==0.0.2',
        'dash-html-components==2.0.0',
        'dash-table==5.0.0',
        'et-xmlfile==1.1.0',
        'fonttools==4.51.0',
        'greenlet==3.0.3',
        'kiwisolver==1.4.5',
        'markupsafe==2.1.5',
        'matplotlib==3.5.3',
        'mordred==1.2.0',
        'networkx==2.8.8',
        'numpy==1.26.4',
        'openpyxl==3.1.2',
        'pandas==2.2.2',
        'pillow==10.3.0',
        #'pycairo==1.23.0',
        'pysocks==1.7.1',
        'python-dateutil==2.9.0.post0',
        'pytz==2024.1',
        'pyyaml==6.0.2rc1',
        'reportlab==3.6.13',
        'scipy==1.13.1',
        'sqlalchemy==2.0.30',
        'tzdata==2024.1',
        'unicodedata2==15.1.0',
        'PyYAML',
        'pandas>=2.0.2',
    ],
    python_requires=">=3.0",
    include_package_data=True,
)
