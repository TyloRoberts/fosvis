import pathlib
from setuptools import find_packages, setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="fosvis",
    version="0.0.14",
    description="A tool for visualizing fosmid DNA sequences.",
    long_description=README,
    long_description_content_type="text/markdown",
    # url="https://github.com/realpython/reader",
    author="Tylo Roberts",
    author_email="tylojroberts@gmail.com",
    license="MIT",
    packages=find_packages(exclude=("tests",)),
    include_package_data=True,
    install_requires=['biopython', 'pandas', 'numpy', 'matplotlib', 'seaborn', 'tqdm']
)
