from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pypaperretriever',
    version='0.1.0',
    author='Joseph Turner',
    author_email='josephisaacturner@gmail.com',
    packages=find_packages(),
    url='http://github.com/josephisaacturner/pypaperretriever',
    license='LICENSE.txt',
    description='An automated tool to download scientific papers.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'pypaperretriever = pypaperretriever.pypaperretriever:main'
        ],
    },
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
