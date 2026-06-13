import glob, os
from PyInstaller.utils.hooks import collect_all, collect_submodules

matplotlib_datas, matplotlib_binaries, matplotlib_hidden = collect_all('matplotlib')
numpy_datas, numpy_binaries, numpy_hidden = collect_all('numpy')
pandas_datas, pandas_binaries, pandas_hidden = collect_all('pandas')

# Manually grab _backend_agg.so
import matplotlib
mpl_backends_dir = os.path.join(os.path.dirname(matplotlib.__file__), 'backends')
agg_binaries = [
    (so, 'matplotlib/backends')
    for so in glob.glob(os.path.join(mpl_backends_dir, '_backend_agg*.so'))
]

a = Analysis(
    ['main.py'],
    pathex=[],
    binaries=(
        matplotlib_binaries + numpy_binaries +
        pandas_binaries  +
        agg_binaries          # <-- explicit injection
    ),
    datas=matplotlib_datas + numpy_datas + pandas_datas,
    hiddenimports=(
        matplotlib_hidden + numpy_hidden + pandas_hidden
        + collect_submodules('matplotlib.backends')
        + collect_submodules('concurrent')
        + collect_submodules('multiprocessing')
        + [
            'matplotlib.backends.backend_agg',
            'matplotlib.backends._backend_agg',
            'concurrent.futures',
            'concurrent.futures._base',
            'concurrent.futures.thread',
            'concurrent.futures.process',
        ]
    ),
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['matplotlib.tests', 'numpy.f2py.tests', 'pandas.tests'],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='RGT',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,          # disable UPX — it corrupts .so files silently
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)