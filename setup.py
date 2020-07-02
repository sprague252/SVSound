import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SVSound3-spraguem", # Replace with your own username
    version="0.0.1",
    author="Mark Sprague",
    author_email="spraguem@ecu.edu",
    description="A package to read Broadcast Wave files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sprague252/SVSound",
    package_dir = {'': 'src'},
    packages=setuptools.find_packages('SVSound', 'SVSound.*'),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    #install_requires=['numpy'],
    python_requires='>=3.6',
)