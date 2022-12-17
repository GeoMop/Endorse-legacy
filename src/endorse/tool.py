import os
from pathlib import Path
import shutil

class workdir:
    """
    Context manager for creation a workspace and change CWD to it.
    """
    def __init__(self, name="sandbox", clean=False):
        self.work_dir = os.path.abspath(name)
        Path(self.work_dir).mkdir(parents=True, exist_ok=True)
        self._clean = clean
        self._orig_dir = os.getcwd()
        os.chdir(self.work_dir)

    def __enter__(self):
        return self.work_dir

    def __exit__(self, type, value, traceback):
        os.chdir(self._orig_dir)
        if self._clean:
            shutil.rmtree(self.work_dir)

