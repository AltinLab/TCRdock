# setup.py
import subprocess
import sys
from setuptools import setup, find_packages
from setuptools.command.install import install


class CustomInstall(install):
    """Run our database installer, then fall back to the normal install."""

    def run(self):
        subprocess.check_call([sys.executable, "download_blast.py"])
        super().run()


setup(
    name="tcrdock",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython==1.79",
        "numpy==1.19.5",
        "pandas==1.3.4",
        "scipy==1.7.0",
        "matplotlib==3.3.4",
    ],
    cmdclass={
        "install": CustomInstall,
    },
)
