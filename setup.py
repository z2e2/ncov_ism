# https://wiki.python.org/moin/Distutils/Tutorial
from setuptools import setup
import ncov_ism

setup(name='ncov_ism',
	version=ncov_ism.__version__,
	description='Informative Subtype Marker (ISM) is an efficient framework for genetic subtyping of a pandemic virus and implement it for SARS-CoV-2, the novel coronavirus that causes COVID-19.',
	url='https://github.com/z2e2/ncov_ism',
	author='Zhengqiao Zhao',
        author_email='zz374@drexel.edu',
	classifiers=['Development Status :: 1 - Beta',
	'Environment :: Console',
	'Intended Audience :: Science/Research',
	'License :: OSI Approved :: BSD License',
	'Natural Language :: English',
	'Operating System :: POSIX :: Linux',
	'Topic :: Scientific/Engineering :: Bio-Informatics'],
        license='BSD',
	packages=['ncov_ism'],
	install_requires=['matplotlib>=2.1.1','pandas>=0.20.3',
                          'numpy>=1.13.3','biopython>=1.71'],
	entry_points = {
        'console_scripts': ['ncov_ism=ncov_ism.ism:main'],
    },
	zip_safe=False)
