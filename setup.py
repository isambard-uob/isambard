"""Setup script for the ISAMBARD."""

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize


def readme():
    """Loads the readme file for AMPAL."""
    with open('README.md', 'r') as inf:
        return inf.read()


setup(name='ISAMBARD',
      version='2.0.1',
      description='A package for biomolecular analysis, modelling and design',
      long_description=readme(),
      long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
      url='https://github.com/isambard-uob/isambard',
      author='Woolfson Group, University of Bristol',
      author_email='chris.wood@bristol.ac.uk',
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      packages=find_packages('src'),
      package_dir={'': 'src'},
      include_package_data=True,
      # This code automatically builds the Cython extensions.
      ext_modules=cythonize(
          [Extension(
              "isambard.specifications.ta_polypeptide",
              [("src/isambard/specifications/ta_polypeptide.pyx")]),
           ]
      ),
      install_requires=[
          'ampal',
          'budeff',
          'Cython',
          'deap',
          'matplotlib',
          'numpy',
      ],
      zip_safe=False,
      )
