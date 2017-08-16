import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except(IOError, ImportError):
    long_description = open('README.md').read()


setuptools.setup(
    name='isambard',
    packages=setuptools.find_packages(),
    # This code automatically builds the Cython extensions.
    ext_modules=cythonize(
        [Extension(
            "isambard_dev.tools.geometry",
            ["isambard_dev/tools/geometry.pyx"]),
         Extension(
             "isambard_dev.buff.calculate_energy",
             ["isambard_dev/buff/calculate_energy.pyx"],
             include_dir=["isambard_dev/buff/"],
             language='c++'),
         Extension(
             "isambard_dev.ampal.specifications.polymer_specs.ta_polypeptide",
             [("isambard_dev/ampal/specifications/polymer_specs/"
               "ta_polypeptide.pyx")]),
         ]),
    include_package_data=True,
    version='2017.2.2',
    description=(
        'ISAMBARD: An open-source computational environment for'
        ' biomolecular analysis, modelling and design'),
    long_description=long_description,
    author='Woolfson Group, University of Bristol',
    author_email='isambardinfo@gmail.com',
    url='https://github.com/woolfson-group/isambard/',
    download_url='https://github.com/woolfson-group/isambard/tarball/2017.2.2',
    keywords=['isambard', 'biomolecule', 'parametric',
              'modelling', 'bristol', 'woolfson'],
    install_requires=[
        'Cython',
        'numpy',
        'requests',
        'SQLAlchemy',
        'networkx',
        'ipywidgets',
        'deap',
        'matplotlib',
        'hypothesis',
        'pyreadline',
        'bs4',
        'parmed',
        'recommonmark',
        'numpydoc',
        'pypandoc'
    ]
)
