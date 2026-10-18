from setup import PlatformWheel
from setuptools.command.bdist_wheel import bdist_wheel


def test_platform_wheel_is_not_marked_pure():
    assert issubclass(PlatformWheel, bdist_wheel)
    assert PlatformWheel.finalize_options is not bdist_wheel.finalize_options
    assert PlatformWheel.get_tag is not bdist_wheel.get_tag
