import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
  name = 'sps',
  version = "1.0",
  author="Daniele Michilli",
  author_email="danielemichilli@gmail.com",
  description="Classifier for single-pulse searches in the fast radio timing field of Astronomy",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/danielemichilli/SpS",
  packages=setuptools.find_packages(),
  classifiers=(
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ),
  install_requires=[
        'python_version<3',
        'cython>=0.28.3',
        'matplotlib>=2.2.2',
        'numpy>=1.14.3',
        'pandas>=0.23.0',
        'scipy>=0.19.1',
        'astropy>=2.0.2',
      ],
)
