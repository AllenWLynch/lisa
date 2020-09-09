

from setuptools import setup

with open('lisa/_version.py', 'r') as v:
    version = v.read().split('=')[-1].strip()

setup(
    name = 'lisa2',
    description = """Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data\n
X. Shirley Liu Lab, 2020""",
    version = '2.0.2',
    url = 'https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6',
    author = 'Allen Lynch',
    author_email = 'alynch@ds.dfci.harvard.edu',
    packages = ['lisa'],
    zip_safe = False,
    scripts = ['bin/lisa'],
    install_requires = [
        'numpy>=1.19,<2',
        'scipy>=2.10.0,<3',
        'h5py>=2.10.0,<3',
        'scikit-learn>=0.23.2,<1'
    ],
    include_package_data=True,
    python_requires='>=3.6,<4',
)