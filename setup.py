from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="dandd",
    version="1.0.0",
    author="Jessica K. Bonnie",
    author_email="jbonnie@example.com",  # Update with actual email
    description="A tool to estimate deltas for sequence sets and answer questions about relative contribution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jessicabonnie/dandd",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy",
        "pandas",
        "biopython",
    ],
    include_package_data=True,
    package_data={
        'dandd': ['lib/*'],
    },
    entry_points={
        'console_scripts': [
            'dandd=dandd.cli:main',
        ],
    },
)