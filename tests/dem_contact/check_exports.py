#!/usr/bin/env python3
"""Check data values, row selection, and every XDMF dataset reference."""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
import math
import h5py
import numpy as np


def verify(prefix):
    expected = {
        'overlap-zero': ('NormalOverlap', [.1]),
        'force-only': ('Normal', [1, 0, 0]),
        'tangential-only': ('Tangential', [0, 1e-15, 0]),
        'dashpot-only': ('DashpotForceTangential', [0, 2, 0]),
    }
    for case, (dataset, values) in expected.items():
        stem = Path(str(prefix) + '-' + case)
        with h5py.File(stem.with_suffix('.h5')) as data:
            assert data['ID1'].size == data['ID2'].size == 1
            np.testing.assert_allclose(data[dataset][:], values, rtol=1e-6, atol=1e-22)
            np.testing.assert_array_equal(data['Rolling'][:], [3, 4, 5])
            assert all(np.isfinite(data[key][:]).all() for key in data)
        for item in ET.parse(stem.with_suffix('.xmf')).findall('.//DataItem'):
            filename, key = item.text.strip().split(':', 1)
            with h5py.File(filename) as data:
                assert key in data, (stem, key)
                assert data[key].size == math.prod(map(int, item.attrib['Dimensions'].split()))
    assert not Path(str(prefix) + '-empty.h5').exists()
    print('CONTACT_EXPORT_PASS cases=5')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    verify(parser.parse_args().prefix)
