from setuptools import setup, find_packages

f = open("README.md", "r")
readme = f.read()
f.close()


setup(
    name="luseesky",
    version="0.0.1",
    # description="",
    long_description=readme,
    author="Christian H. Bye",
    author_email="chbye@berkeley.edu",
    url="https://github.com/christianhbye/lusee_sky_simulations",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "pyuvsim",
        "pyuvdata",
        "jupyter",
        "ipykernel",
    ],
    extras_require={"tests": "pytest", "style": ["black", "flake8"]},
    include_package_data=True,
    python_requires=">=3.7",
    # license="MIT",
    # classifiers=[
    #    "Intended Audience :: End Users/Desktop",
    #    "License :: OSI Approved :: MIT License",
    #    "Programming Language :: Python :: 3",
    #    "Topic :: Utilities",
    # ],
)
