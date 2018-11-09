#!/usr/bin/env python
"""
Build Python interface for veaslc
"""
import os
import sys
import subprocess
import re

from setuptools import setup, find_packages

def git_version():
    # Full version includes the Git commit hash
    full_version = subprocess.check_output('git describe --dirty', shell=True).decode("utf-8").strip(" \n")

    # Python standardized version in form major.minor.patch.dev<build>
    version_regex = re.compile(r"v?(\d+\.\d+\.\d+(-\d+)?).*")
    match = version_regex.match(full_version)
    if match:
        std_version = match.group(1).replace("-", ".dev")
    else:
        raise RuntimeError("Failed to parse version string %s" % full_version)

    return full_version, std_version

def git_timestamp():
    return subprocess.check_output('git log -1 --format=%cd', shell=True).decode("utf-8").strip(" \n")

def set_metadata(module_dir, version_str, timestamp_str):
    vfile = open(os.path.join(module_dir, "oxasl_ve", "_version.py"), "w")
    vfile.write("__version__ = '%s'\n" % version_str)
    vfile.write("__timestamp__ = '%s'\n" % timestamp_str)
    vfile.close()

# Read in requirements from the requirements.txt file.
with open('requirements.txt', 'rt') as f:
    requirements = [l.strip() for l in f.readlines()]

# Generate a list of all of the packages that are in your project.
packages = find_packages()

rootdir = os.path.join(os.path.abspath(os.path.dirname(__file__)))
_, stdv = git_version()
timestamp = git_timestamp()
set_metadata(rootdir, stdv, timestamp)

desc = "Python interface for vessel-encoded ASL"
longdesc = desc

# setup parameters
setup(name='oxasl_ve',
      version=stdv,
      description=desc,
      long_description=longdesc,
      author='Michael Chappell, Martin Craig',
      author_email='martin.craig@eng.ox.ac.uk',
      packages=packages,
      include_package_data=True,
      data_files=[],
      setup_requires=[],
      install_requires=[],
      entry_points={
      },
      classifiers=["Programming Language :: Python :: 2.7",
                   "Development Status:: 3 - Alpha",
                   'Programming Language :: Python',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: Microsoft :: Windows',
                   'Operating System :: POSIX',
                   "Intended Audience :: Education",
                   "Intended Audience :: Science/Research",
                   "Intended Audience :: End Users/Desktop",
                   "Topic :: Scientific/Engineering :: Bio-Informatics",])
