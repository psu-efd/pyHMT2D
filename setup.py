#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Modified from https://raw.githubusercontent.com/navdeep-G/setup.py/master/setup.py

# Note: To use the 'upload' functionality of this file, you must:
#   $ pip install build twine

import io
import os
import subprocess
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'pyHMT2D'
DESCRIPTION = 'Two-dimensional hydraulic modeling tools in Python'
URL = 'https://github.com/psu-efd/pyHMT2D'
EMAIL = 'xiaofengliu19@gmail.com'
AUTHOR = 'Xiaofeng Liu'
REQUIRES_PYTHON = '>=3.8.0, <3.14.0'
VERSION = '2.0.0'

# What packages are required for this module to be executed?
REQUIRED = [
    'numpy', 'vtk>=9.0.2', 'h5py', 'scipy', 'affine', 'meshio', 'scikit-optimize',
    'rasterio', 'pywin32'
]

# What packages are optional?
EXTRAS = {
    'ai': [
        'mcp>=1.0.0',           # Model Context Protocol server (FastMCP)
    ],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


class UploadCommand(Command):
    """Support setup.py upload.

    Usage:
        python setup.py upload          # build and upload to PyPI
        python setup.py upload --test   # build and upload to Test PyPI
    """

    description = 'Build and publish the package using build + twine.'
    user_options = [
        ('test', 't', 'Upload to Test PyPI instead of production PyPI'),
    ]

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        self.test = False

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
            rmtree(os.path.join(here, 'build'))
        except OSError:
            pass

        self.status('Building source and wheel distribution…')
        subprocess.check_call([sys.executable, '-m', 'build'])

        if self.test:
            self.status('Uploading to Test PyPI via Twine…')
            subprocess.check_call([
                sys.executable, '-m', 'twine', 'upload',
                '--repository', 'testpypi', 'dist/*',
            ])
        else:
            self.status('Uploading to PyPI via Twine…')
            subprocess.check_call([
                sys.executable, '-m', 'twine', 'upload', 'dist/*',
            ])

        self.status('Pushing git tags…')
        tag = 'v{0}'.format(about['__version__'])
        subprocess.check_call(['git', 'tag', tag])
        subprocess.check_call(['git', 'push', '--tags'])

        sys.exit()


# Where the magic happens:
setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    # If your package is a single module, use this instead of 'packages':
    # py_modules=['mypackage'],

    entry_points={
        'console_scripts': ['hmt-mcp-server=pyHMT2D.AI_Tools.mcp_server:main',
                            'hmt-cli=pyHMT2D.AI_Tools.cli_runner:main',
                            ],
    },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license='MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering',
        'Topic :: Utilities'
    ],
    # $ setup.py publish support.
    cmdclass={
        'upload': UploadCommand,
    },
)
