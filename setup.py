from setuptools import setup

setup(name='skmer',
      version='1.0.0',
      description='Assembly-free and alignment-free tool for estimating genomic distances between genome-skims',
      author='Shahab Sarmashghi',
      author_email='ssarmash@ucsd.edu',
      license='BSD-3-Clause',
      url='https://github.com/shahab-sarmashghi/Skmer',
      packages=['skmer'],
      package_dir={'skmer': 'skmer'},
      install_requires=['numpy>=1.13.1', 'scipy>=0.19.1', 'pandas>=0.20.3'],
      provides=["skmer"],
      entry_points={
            'console_scripts': ['skmer=skmer.__main__:main']
      },
      classifiers=["Environment :: Console",
                   "Intended Audience :: Developers",
                   "Intended Audience :: Science/Research",
                   "Operating System :: Unix",
                   "Programming Language :: Python",
                   "Programming Language :: Python :: 2",
                   "Programming Language :: Python :: 3",
                   "Topic :: Scientific/Engineering :: Bio-Informatics"],
      )

