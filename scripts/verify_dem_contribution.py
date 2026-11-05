#!/usr/bin/env python3
"""Build and verify the DEM contribution in a fresh, retained output directory."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sources(root, output):
    result = subprocess.run(['git', 'ls-files', '-z', '--cached', '--others',
                             '--exclude-standard'], cwd=root, check=True,
                            capture_output=True)
    paths = sorted(set(result.stdout.decode().rstrip('\0').split('\0')))
    return {name: (hashlib.sha256(str((root/name).readlink()).encode()).hexdigest()
                   if (root/name).is_symlink() else digest(root/name))
            for name in paths if name and output not in (root/name).parents
            and ((root/name).is_symlink() or (root/name).is_file())}


def replace_parameter(text, name, value):
    lines = text.splitlines()
    matches = [i for i, line in enumerate(lines)
               if re.search(r'\b'+re.escape(name)+r':', line)]
    if len(matches) != 1:
        raise RuntimeError('Expected exactly one input parameter: '+name)
    i = matches[0]
    lines[i] = str(value)+' '+lines[i].split(None, 1)[1]
    return '\n'.join(lines)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--full-examples', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build, out = args.build.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    summary = {'classification': 'FAIL', 'checks': [], 'source': str(root),
               'build': str(build), 'full_examples': args.full_examples,
               'timing': 'NOT_A_PERFORMANCE_MEASUREMENT'}

    def save():
        (out/'SUMMARY.json').write_text(json.dumps(summary, indent=2)+'\n')

    def run(name, command, cwd=out, expected=None, timeout=120):
        row = {'name': name, 'command': list(map(str, command)), 'cwd': str(cwd)}
        summary['checks'].append(row)
        save()
        try:
            result = subprocess.run(row['command'], cwd=cwd, capture_output=True,
                                    text=True, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            row['error'] = 'timeout'
            (out/(name+'.log')).write_bytes((exc.stdout or b'')+(exc.stderr or b''))
            save()
            raise RuntimeError(name+' timed out') from exc
        log = result.stdout+result.stderr
        (out/(name+'.log')).write_text(log)
        row['returncode'] = result.returncode
        save()
        if result.returncode:
            raise RuntimeError(name+' failed; see '+str(out/(name+'.log')))
        if expected and expected not in log:
            raise RuntimeError(name+' did not report '+expected)
        print(name+': PASS', flush=True)
        return log

    try:
        cache_file = build/'CMakeCache.txt'
        if not cache_file.is_file():
            raise RuntimeError('Configure the build before verification: '+str(build))
        cache = dict(line.split('=', 1) for line in cache_file.read_text().splitlines()
                     if '=' in line and not line.startswith(('#', '//')))
        if Path(cache.get('CMAKE_HOME_DIRECTORY:INTERNAL', '')).resolve() != root:
            raise RuntimeError('Build was configured for a different source checkout')
        for option in ['A_WITH_TDEM', 'A_WITH_DEM_CONTACT_REGRESSIONS', 'A_WITH_EXPERIMENTAL_VV_COUPLING']:
            if cache.get(option+':BOOL') != 'ON':
                raise RuntimeError('Configure '+option+'=ON')
        if args.jobs < 1:
            raise RuntimeError('--jobs must be positive')
        summary['source_head'] = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip()
        summary['source_status'] = subprocess.check_output(
            ['git', 'status', '--porcelain=v1', '-uall'], cwd=root, text=True)
        summary['source_sha256'] = sources(root, out)
        run('configure', ['cmake', '-S', root, '-B', build], timeout=120)
        run('build', ['cmake', '--build', build, '--target',
                      'dem_contribution_verification_binaries',
                      '--parallel', args.jobs], timeout=1800)
        listing = run('ctest-list', ['ctest', '--test-dir', build, '-R', '^dem_',
                                    '--show-only=json-v1'])
        expected_tests = {'dem_contact_cpu', 'dem_contact_gpu', 'dem_vv_cpu', 'dem_vv_gpu',
                          'dem_legacy_coupling_cpu', 'dem_legacy_coupling_gpu',
                          'dem_coupling_smoke', 'dem_test_force_api_cpu',
                          'dem_test_periodic_pid_cpu', 'dem_test_periodic_pid_gpu',
                          'dem_test_shearbox_gpu'}
        actual = {test['name'] for test in json.loads(listing)['tests']}
        if not expected_tests <= actual:
            raise RuntimeError('Missing CTest registrations: '+str(expected_tests-actual))
        summary['ctest_names'] = sorted(actual)
        run('ctest', ['ctest', '--test-dir', build, '-R', '^dem_',
                      '--no-tests=error', '--output-on-failure'], timeout=300)
        tests = build/'tests/dem_contact'
        for backend in ['cpu', 'gpu']:
            run('exports-'+backend, [sys.executable, root/'tests/dem_contact/check_exports.py',
                                     tests/('contact_'+backend)], expected='CONTACT_EXPORT_PASS')
        run('histories', [sys.executable, root/'tests/dem_contact/check_histories.py',
                          tests, out/'histories'], expected='HISTORY_REGRESSION_PASS', timeout=300)
        case = out/'saveload'; case.mkdir()
        run('saveload', [build/'tdem/test_saveload'], cwd=case,
            expected='PASS - the reloaded domain is bit-for-bit the same')
        summary['inputs'] = {}
        for integrator in [0, 1]:
            # test_cu_02 runs long enough for its analytic turntable orbit to be
            # discriminating: 0.1 s separates a correct run from one whose
            # rotational path is disabled by nearly three orders of magnitude.
            for target, extra in [('test_cu_01', ['1', '.05']),
                                  ('test_cu_02', ['1', '.1']),
                                  ('test_cu_03', ['1', '.05', '128', '.01']),
                                  ('test_cu_04', ['case', '1', '4.0', '0.0'])]:
                name = target+'-'+str(integrator)
                case = out/name; case.mkdir()
                if target == 'test_cu_04':
                    shutil.copyfile(root/'tests/dem_contact/inputs/test_cu_04.inp', case/'case.inp')
                    summary['inputs'][name] = digest(case/'case.inp')
                run(name, [build/'tdem'/target, *extra, str(integrator)], cwd=case,
                    expected='DEM_BENCHMARK_INTEGRATOR velocity_verlet='+str(integrator), timeout=600)
                if target == 'test_cu_02':
                    # Two spinning bodies and a closed-form orbit: this is the
                    # analytic check of the rotational integration path.
                    run('turntable-'+str(integrator),
                        [sys.executable, root/'tests/dem_contact/check_turntable.py',
                         case/'test_cu_02_trayectory.res'],
                        expected='TURNTABLE_ROLLING_PASS')
        for target, input_name in [('PID_v2', 'test_pid_v2.inp'),
                                   ('TX_PID_periodic', 'test_tx_pid_periodic.inp'),
                                   ('TX_PID_periodic_gpu', 'test_tx_pid_periodic.inp'),
                                   ('shearbox_LE_const_vol', 'test_shearbox_LE.inp')]:
            case = out/target; case.mkdir()
            text = (root/'tdem'/input_name).read_text()
            if not args.full_examples:
                if target == 'shearbox_LE_const_vol':
                    # Retain both phase durations and dt; reduce only the box.
                    for parameter, value in [('Lx', .4), ('Ly', .4), ('Lz', 1.4)]:
                        text = replace_parameter(text, parameter, value)
                else:
                    for parameter, value in [('isotropicEnd', .005), ('shearEnd', .01)]:
                        text = replace_parameter(text, parameter, value)
            (case/'case.inp').write_text(text)
            summary['inputs'][target] = digest(case/'case.inp')
            message = {'shearbox_LE_const_vol': 'All phases completed successfully.',
                       'TX_PID_periodic_gpu': 'TX_PID_periodic (GPU) finished successfully.'}.get(
                           target, target+' finished successfully.')
            run(target, [build/'tdem'/target, 'case', '1'], cwd=case,
                expected=message, timeout=900)
        # Check saved numerical artifacts, not just banners or process completion.
        import h5py
        import numpy as np
        files = 0
        for path in out.rglob('*'):
            if path.suffix in {'.h5', '.hdf5'}:
                with h5py.File(path) as data:
                    def finite(name, item):
                        if isinstance(item, h5py.Dataset) and item.dtype.kind in 'fc':
                            if not np.isfinite(item[()]).all():
                                raise RuntimeError('Nonfinite data in '+str(path)+':'+name)
                    data.visititems(finite)
                files += 1
            elif path.suffix == '.res':
                for token in path.read_text().split():
                    try:
                        number = float(token)
                    except ValueError:
                        continue
                    if not math.isfinite(number):
                        raise RuntimeError('Nonfinite result in '+str(path))
        if files == 0:
            raise RuntimeError('No numerical artifacts produced')
        summary['finite_hdf5_files'] = files
        summary['binaries'] = {str(path.relative_to(build)): digest(path)
                               for directory in [tests, build/'tdem']
                               for path in directory.iterdir()
                               if path.is_file() and path.stat().st_mode & 0o111
                               and path.suffix == ''}
        if summary['source_sha256'] != sources(root, out):
            raise RuntimeError('Source files changed during verification')
        summary['classification'] = 'PASS_BOUNDED_REGRESSIONS'
        save()
        print('DEM_CONTRIBUTION_VERIFICATION_PASS '+str(out/'SUMMARY.json'))
        return 0
    except Exception as exc:
        summary['error'] = str(exc)
        save()
        print('DEM_CONTRIBUTION_VERIFICATION_FAIL: '+str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
