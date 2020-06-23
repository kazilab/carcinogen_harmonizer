"""
Setup kFDR
"""

import re
import sys
try:
    from setuptools import setup
except ImportError:
    raise ImportError("Please install `setuptools`")

if not (sys.version_info[0] == 3 and sys.version_info[1] >= 7):
    raise RuntimeError(
                'kFDR requires Python 3.7 or higher.\n'
                'You are using Python {0}.{1}'.format(
                    sys.version_info[0], sys.version_info[1])
                )

packagedata={'kFDR': ['data/*']}
# main setup command
setup(
    name='kFDR',
    version='0.0.1',
    author='Julhash Kazi',
    author_email='jk@kazilab.se',
    url='https://www.kazilab.se',
    description='To calculate p-value, fold changes and adjusted p-value',
    license='GPL',
    install_requires=[
        'statsmodels>=0.11.1',
        'pandas>=0.25.3',
        'numpy>=1.18.2',
        'scipy>=1.4.1',
        ],
    platforms='Mac OS X',
    packages=['kFDR'],
    package_dir={'kFDR': 'src/kFDR'},
    package_data=packagedata
)
