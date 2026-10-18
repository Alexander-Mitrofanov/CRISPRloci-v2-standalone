from setuptools import setup, find_packages
setup(name="phigaro-v2", version="2.4.1", packages=find_packages(),
      description="Bounded offline Phigaro with source-bound evidence",
      license="MIT", python_requires=">=3.11",
      package_data={"phigaro":["*.json"]},
      install_requires=["numpy==2.2.6", "biopython==1.85", "six==1.17.0"],
      entry_points={"console_scripts":["phigaro-v2=phigaro.bounded:main"]})
