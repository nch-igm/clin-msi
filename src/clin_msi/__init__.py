import os

version_file = os.path.join(os.path.dirname(__file__), "_version.py")

try:
    from ._version import __version__
except ModuleNotFoundError:
    raise RuntimeError(
        f'There was no file {version_file} with given version found. '
        f'Run in command line: echo -n __version__ = \\"$(dunamai from git --style semver)\\" > {version_file}'
    )
