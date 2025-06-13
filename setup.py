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
    packages=["tcrdock"],
    install_requires=[
        "biopython",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
    ],
    cmdclass={
        "install": CustomInstall,
    },
    package_data={
        "tcrdock": [
            "**/*.alfas",
            "**/*.fasta",
            "**/*.txt",
            "**/*.tsv",
            "**/tcr/*.pdb",
            "**/tcr/*.json",
            "**/ternary/*.pdb",
            "**/ternary/*.json",
            "**/pmhc/*.pdb",
            "**/pmhc/*.json",
        ]
    },
)
