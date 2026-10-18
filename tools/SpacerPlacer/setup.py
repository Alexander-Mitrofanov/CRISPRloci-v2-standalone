"""Compatibility shim for legacy installers; metadata lives in pyproject.toml."""

from setuptools import setup
from setuptools.command.bdist_wheel import bdist_wheel


class PlatformWheel(bdist_wheel):
    """Tag wheels by platform because bundled MAFFT converters are native binaries."""

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        _, _, platform_tag = super().get_tag()
        return 'py3', 'none', platform_tag


if __name__ == "__main__":
    setup(cmdclass={'bdist_wheel': PlatformWheel})
