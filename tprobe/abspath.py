import os
from pathlib import Path, _posix_flavour, _windows_flavour


class AbsPath(Path):
    _flavour = _windows_flavour if os.name == 'nt' else _posix_flavour

    @property
    def abspath(self):
        return self.resolve().__str__()

    @property
    def str(self):
        return self.abspath

