import os
from setuptools import setup, find_packages

setup(name='pyatactools',
      version='0.0.1',
      #packages=find_packages(),
      description='pyatactools is a Python module to analyze ATAC-seq NGS data',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=['pyatactools'],
      package_data={"pyatactools":['data/*']},
      scripts=['scripts/pyatac_preprocess.py', 'scripts/pyatac_align.py', 'scripts/pyatac_profiles.py'],
      install_requires=['pysam', 'pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Environment :: Console',
      ],
      long_description="""

pyatactools is a Python module to analyze ATAC-seq NGS data

 Contact
=============

If you have any questions or comments about pychiptools, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
