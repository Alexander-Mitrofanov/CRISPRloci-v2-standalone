from setuptools import setup, Extension

# Compiled as a setuptools extension to retain wheel/platform identities, loaded
# only by ctypes (it is a C ABI library, not an importable Python extension).
setup(ext_modules=[Extension('isescan_v2.native.ssw201507.libssw',
    sources=['isescan_v2/native/ssw201507/ssw.c'], extra_compile_args=['-O3'])])
