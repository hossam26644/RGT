'''Golden-output regression test: RGT on the first 750 demo reads with Demo/HD.json (~1 min).

Run:        python -m pytest tests            (from the repo root)
Regenerate: python main.py -i tests/data/demo750 -o tests/golden/demo_HD -s Demo/HD.json
'''
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
from matplotlib.image import imread

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests/golden/demo_HD"


def files(directory):
    return sorted(p.relative_to(directory) for p in directory.rglob("*") if p.is_file())


def xlsx_parts(path):
    # docProps/core.xml holds creation timestamps; every other part (values, fills, rich text) must match
    with zipfile.ZipFile(path) as archive:
        return {name: archive.read(name) for name in archive.namelist() if name != "docProps/core.xml"}


def assert_same_outputs(output, golden):
    assert files(output) == files(golden)
    for name in files(golden):
        new, old = output / name, golden / name
        if name.suffix == ".xlsx":
            assert xlsx_parts(new) == xlsx_parts(old), name
        elif name.suffix == ".png":
            assert np.array_equal(imread(new), imread(old)), name
        else:
            assert new.read_bytes() == old.read_bytes(), name


def test_demo_hd(tmp_path):
    output = tmp_path / "out"
    subprocess.run([sys.executable, "main.py", "-i", "tests/data/demo750", "-o", str(output),
                    "-s", "Demo/HD.json"], cwd=ROOT, check=True)
    assert_same_outputs(output, GOLDEN)
