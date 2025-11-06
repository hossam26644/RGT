import json
from pathlib import Path
import pytest

@pytest.fixture(scope='session')
def repo_root():
    return Path(__file__).resolve().parents[1]

def load_settings(path):
    with open(path) as fh:
        return json.load(fh)

@pytest.fixture(scope='session')
def sample_settings(repo_root):
    return load_settings(repo_root / 'tests' / 'data' / 'sample_settings.json')

@pytest.fixture(scope='session')
def sample_fastq(repo_root):
    return repo_root / 'tests' / 'data' / 'sample.fastq.txt'

@pytest.fixture(scope='session')
def invalid_fastq(repo_root):
    return repo_root / 'tests' / 'data' / 'invalid.fastq'

@pytest.fixture
def temp_output_dir(repo_root, tmp_path):
    output_dir = tmp_path / 'test_output'
    yield output_dir
    # tmp_path handles cleanup automatically

@pytest.fixture
def no_plot(monkeypatch):
    try:
        from graphsplotter import plotter as plotter_module
        monkeypatch.setattr(plotter_module, 'plot_graphs', lambda *a, **k: None)
    except (ImportError, AttributeError):
        pytest.skip("Skipping plotter mock due to missing module or attribute")
    try:
        from graphsplotter import table_2d_plotter as t2
        monkeypatch.setattr(t2, 'plot_2d_table', lambda *a, **k: None)
    except (ImportError, AttributeError):
        pytest.skip("Skipping 2D plotter mock due to missing module or attribute")
    try:
        from graphsplotter import plot_3D
        monkeypatch.setattr(plot_3D, 'plot_3D', lambda *a, **k: None)
    except (ImportError, AttributeError):
        pytest.skip("Skipping 3D plotter mock due to missing module or attribute")
    # patch the name imported directly into rgt module so rgt.rgt() uses the mocked function
    try:
        import rgt
        monkeypatch.setattr(rgt, 'plot_graphs', lambda *a, **k: None)
    except Exception:
        # non-fatal: if rgt cannot be imported in the test environment, ignore
        pass

@pytest.fixture
def no_excel_save(monkeypatch):
    try:
        from excelexporter.ExcelExport import ExcelWriter
        monkeypatch.setattr(ExcelWriter, 'save_file', lambda self, *a, **k: None)
    except (ImportError, AttributeError):
        pytest.skip("Skipping Excel save mock due to missing module or attribute")