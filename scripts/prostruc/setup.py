from setuptools import setup, find_packages
import os

def readme():
    with open('README.md', encoding='utf-8') as f:
        return f.read()

setup(
    name='prostruc',
    version='0.0.1',
    author='The Prostruc Team',
    author_email='senawilson123@gmail.com',
    description='Prostruc: A Comprehensive Command-line Tool for Protein Structure Prediction and Validation',
    long_description=readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/omicscodeathon/prostruct',
    packages=find_packages(),
    install_requires=[
        'biopython>=1.76',
        'docker>=3.0,<5.0',
        'reportlab>=3.5,<4.0',
    ],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
    ],
    entry_points={
        'console_scripts': [
            'prostruc=prostruc.app:main',
        ],
    },
    include_package_data=True,
    zip_safe=False,
    license='MIT',
    python_requires='>=3.6',
)
