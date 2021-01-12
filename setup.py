

from setuptools import setup, find_packages

setup(
    name = 'lisa2',
    description = """Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data\n
X. Shirley Liu Lab, 2020""",
    version = '2.2.0',
    url = 'https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6',
    author = 'Allen Lynch',
    author_email = 'alynch@ds.dfci.harvard.edu',
    packages = find_packages(),
    zip_safe = True,
    scripts = ['bin/lisa'],
    install_requires = [
        'numpy>=1.19,<2',
        'scipy>=1.5,<2',
        'h5py>=2.10.0,<3',
        'scikit-learn>=0.23.2,<1',
        'pyBigWig>=0.3.17,<1'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    python_requires='>=3.6,<4',
)