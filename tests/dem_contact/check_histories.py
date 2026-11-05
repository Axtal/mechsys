#!/usr/bin/env python3
"""Run bounded CPU/CUDA history, geometry, and checkpoint comparisons."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess


def compare(a, b, path='root', stats=None):
    if stats is None:
        stats = {'numbers': 0, 'max_absolute_difference': 0.0}
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path+'.'+key, stats)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+f'[{i}]', stats)
    elif path.rsplit('.', 1)[-1] in {'step', 'id', 'tag', 'a', 'b', 'contact', 'family', 'key', 'Nc'}:
        assert type(a) is type(b) is int and a == b, (path, a, b)
    else:
        assert math.isfinite(a) and math.isfinite(b), path
        error = abs(a-b)
        assert error <= 1e-10 + 1e-9*max(abs(a), abs(b)), (path, a, b, error)
        stats['numbers'] += 1
        stats['max_absolute_difference'] = max(stats['max_absolute_difference'], error)
    return stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary_directory', type=Path)
    p.add_argument('fresh_output', type=Path)
    args = p.parse_args()
    out = args.fresh_output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binaries = args.binary_directory.resolve()
    results = {}

    def run(label, backend, mode, scenario, law, steps, checkpoint):
        command = [str(binaries/f'dem_test_history_{backend}'), mode,
                   str(scenario), str(law), '1', str(steps), str(checkpoint)]
        env = dict(os.environ, MECHSYS_V2_FEATURE_ORACLE='1')
        result = subprocess.run(command, cwd=out, env=env, capture_output=True,
                                text=True, timeout=120)
        (out/(label+'.log')).write_text(result.stdout+result.stderr)
        assert result.returncode == 0, (label, result.returncode)
        rows = [json.loads(line.split(' ', 1)[1]) for line in result.stdout.splitlines()
                if line.startswith('V2_STATE ')]
        assert rows, label
        # Final callback repeats the last physical boundary only if output cadence
        # already sampled it. Compare by actual integration step in either case.
        return {r['step']: r for r in rows}

    for scenario, law, steps in [(0, 0, 32), (0, 1, 32), (5, 0, 32),
                                 (11, 0, 8), (12, 0, 8), (13, 0, 8), (14, 0, 8)]:
        label = f'scenario{scenario}-law{law}'
        cpu = run(label+'-cpu', 'cpu', 'full', scenario, law, steps, out/'unused')
        gpu = run(label+'-gpu', 'gpu', 'full', scenario, law, steps, out/'unused')
        results[label] = compare(list(cpu.values()), list(gpu.values()))
    for backend in ['cpu', 'gpu']:
        label = 'checkpoint-'+backend
        checkpoint = out/(label+'.json')
        full = run(label+'-full', backend, 'full', 0, 0, 32, checkpoint)
        run(label+'-split', backend, 'split', 0, 0, 16, checkpoint)
        continued = run(label+'-continued', backend, 'continue', 0, 0, 32, checkpoint)
        results[label] = compare([full[k] for k in continued], list(continued.values()))
    (out/'SUMMARY.json').write_text(json.dumps(results, indent=2)+'\n')
    print('HISTORY_REGRESSION_PASS comparisons='+str(len(results)))


if __name__ == '__main__':
    main()
