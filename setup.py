import setuptools
from setuptools import find_packages
__version__="0.1.0"

with open("README.md", "r") as fh:
    long_description = fh.read()


# def packages_list():
#     packages = []
#     packages.extend(find_packages(where='src/endorse'))
#     print("Detected packages: ", packages)
#     [*find_packages(), 'endorse.mesh', 'endorse.common'],
#     return packages

setuptools.setup(
    name="endorse1",
    version=__version__,
    license='GPL 3.0',
    description='Stochastic model of excavation damage zone and its safety indicators.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Jan Brezina',
    author_email='jan.brezina@tul.cz',
    url='https://github.com/flow123d/swrap',
    download_url='https://github.com/flow123d/swrap/archive/v{__version__}.tar.gz',
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Scientific/Engineering',
    ],

    keywords=[
        'excavation damage zone', 'stochastic simulation', 'radioactive waste repository',
    ],
    # include_package_data=True, # package includes all files of the package directory
    zip_safe=False,
    install_requires=['numpy>=1.13.4', 'pandas', 'scipy', 'bih', 'gmsh>=4.10.4', 'pyyaml',
                      'pyyaml-include', 'matplotlib', 'pyvista'],
    python_requires='>=3',


    packages=['endorse', 'endorse.common', 'endorse.mesh', 'endorse.mlmc'],
    package_dir={
        'endorse': 'src/endorse',
    },
    package_data={
        'endorse/flow123d_inputs': ['*.yaml']
    },
    entry_points={
        'console_scripts': ['endorse_gui=endorse.gui.app:main']
    }
)



        

        
