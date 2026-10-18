"""Build the dedicated offline context distribution.

Modified for CRISPRloci 2026-10-12. Upstream GPL-3.0 source remains in this
repository; only the bounded context package is installed by this distribution.
"""
from pathlib import Path
from setuptools import setup

ROOT = Path(__file__).resolve().parent
requirements = [line.strip() for line in (ROOT / 'requirements-context.lock').read_text().splitlines()
                if line.strip() and not line.lstrip().startswith('#')]
setup(
    name='crisprloci-defense-context',
    version='2.0.1',
    description='Offline source-bound DefenseFinder HMM context',
    python_requires='>=3.12',
    license='GPL-3.0',
    license_files=['LICENSE', 'third_party/*LICENSE'],
    packages=['defense_context'],
    install_requires=requirements,
    entry_points={'console_scripts': ['defense-context=defense_context.pipeline:main']},
    package_data={'defense_context': ['models.lock.json']},
)
