from setuptools import setup
from setuptools import find_packages

with open("README.md","r", encoding="utf-8") as fh:
    long_description = fh.read()
    
setup(
    name="sedcreator",
    version='0.8.6',
    description='sedcreator is a package that has two main classes, SedFluxer and SedFitter. SedFluxer performs aperture photometry on a given image, coordinates and aperture size. It has a number of functions to print useful information and to plot the image together with apertures. SedFitter fits observations to a grid of models following the Zhang and Tan (2018) radiative transfer models.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        'Topic :: Scientific/Engineering :: Astronomy',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],
    py_modules=['sedcreator'],
    url='https://github.com/fedriani/sedcreator',
    project_urls={"Bug Tracker":"https://github.com/fedriani/sedcreator/issues"},
    author='Ruben Fedriani',
    author_email='ruben.fedriani@gmail.com',
    license='GNU General Public License v3.0',
    install_requires=[
        'astropy>=4.3.1',
        'astroquery>=0.4.3',
        'matplotlib',
        'numpy',
        'photutils>=1.3.0',
        'scipy>=1.5.0',
        'setuptools',
        'tqdm'
        ],
    packages=find_packages(),
    python_requires=">=3.6",
    include_package_data=True,
    package_data={'sedcreator': ['Model_SEDs/*']},
)
