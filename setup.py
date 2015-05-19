from setuptools import setup

from royfit import version

setup(name='royfit',
      version=version,
      url='http://github.com/tamasgal/royfit/',
      description='Restless Oyster Fit - Fast muon reconstruction for KM3NeT',
      author='Tamas Gal',
      author_email='tgal@km3net.de',
      packages=['royfit'],
      include_package_data=True,
      platforms='any',
      install_requires=[
          'numpy',
          'km3pipe',
          'iminuit',
          'matplotlib',
      ],
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
      ],
)

__author__ = 'Tamas Gal'
