# Usage: python angles.py [run_id]

import os
import datetime
import itertools

import numpy
import pandas

import hxform
import utilrsw

log = utilrsw.logger(console_format='%(message)s')


def libs(frame_in, frame_out, excludes=None):
  lib_infos = hxform.lib_info()
  libs_avail = []
  for lib, lib_info in lib_infos.items():
    if excludes is not None and lib in excludes:
      continue
    if frame_in not in lib_info['frames'] or frame_out not in lib_info['frames']:
      wmsg = f"{frame_in} and {frame_out} are both not available in {lib}"
      wmsg += ", skipping this transform pair comparison."
      log.warning(wmsg)
      continue
    libs_avail.append(lib)
  return libs_avail


def angles(to, tf, axis, delta, libs, hxform_args):

  if axis not in ['x', 'y', 'z']:
    raise ValueError(f"Invalid axis: {axis}. Must be one of 'x', 'y', or 'z'.")

  t = hxform.timelib.ints_list(to, tf, delta)
  t = numpy.array(t, dtype=numpy.int32)
  t = t[:, 0:6]
  t_dts = hxform.timelib.ints2datetime(t)
  Δθ = numpy.full((t.shape[0], len(libs)), numpy.nan)
  δψ = numpy.full((t.shape[0], len(libs)), numpy.nan) # Uncertainty in ψ ≔ Δθ

  for i, lib in enumerate(libs):
    #log.info(f"Processing {lib}...")
    hxform_args['lib'] = lib

    if axis == 'x':
      p_in = numpy.array([1., 0., 1.])
    if axis == 'y':
      p_in = numpy.array([0., 1., 0.])
    if axis == 'z':
      p_in = numpy.array([0., 0., 1.])

    p_out = hxform.transform(p_in, t, **hxform_args)

    n = numpy.dot(p_out, p_in)
    d = numpy.linalg.norm(p_out, axis=1)*numpy.linalg.norm(p_in)
    Δθ[:, i] = (180.0/numpy.pi)*numpy.arccos(n/d)

    ε = numpy.spacing(numpy.abs(Δθ[:, i]).max())
    if lib == 'sscweb':
      # values reported to two decimal places => max uncertainty of 0.005
      ε = 0.005
    δψ[:, i] = _angle_uncert(p_in, p_out, n, ε)

    if lib.startswith('spiceypy1'):
      years = numpy.array([dt.year for dt in t_dts])
      mask = (years < 1990) | (years >= 2020)
      Δθ[mask, i] = numpy.nan

    if lib.startswith('spiceypy2'):
      years = numpy.array([dt.year for dt in t_dts])
      mask = (years < 1990) | (years >= 2030)
      Δθ[mask, i] = numpy.nan

  return t_dts, Δθ, δψ


def diffs(dfs, xform, index, libs_avail):

  # Compute diffs DataFrame
  columns_diff = libs_avail.copy()
  columns_diff.remove('geopack_08_dp')
  df = pandas.DataFrame(numpy.nan, index=index, columns=columns_diff)
  if 'geopack_08_dp' in libs_avail:
    for lib in libs_avail:
      if lib != 'geopack_08_dp':
        diff = dfs[xform]['values'][lib] - dfs[xform]['values']['geopack_08_dp']
        df[lib] = diff

  max_min = dfs[xform]['values'].max(axis=1) - dfs[xform]['values'].min(axis=1)
  df['|max-min|'] = numpy.abs(max_min)

  return df


def _angle_uncert(p_in, p_out, n, ε):
  """Estimate uncertainty in angular difference Δθ using propagation of error.
  ψ ≔ Δθ = acos(n/d)
  assume d = 1

  dψ = -dn/(sqrt(1-n^2))

  n = dot(p_out, p_in) = p_out_x*p_in_x + p_out_y*p_in_y + p_out_z*p_in_z

  dn = ( |∂n/∂p_out_x| * Δp_out_x
       + |∂n/∂p_out_y| * Δp_out_y
       + |∂n/∂p_out_z| * Δp_out_z
       + |∂n/∂p_in_x| * Δp_in_x
       + |∂n/∂p_in_y| * Δp_in_y
       + |∂n/∂p_in_z| * Δp_in_z)

  ∂n/∂p_out_x = p_in_x
  ∂n/∂p_out_y = p_in_y
  ∂n/∂p_out_z = p_in_z
  ∂n/∂p_in_x  = p_out_x
  ∂n/∂p_in_y  = p_out_y
  ∂n/∂p_in_z  = p_out_z

  Let ε ≔ Δp_out_i = Δp_in_i for i = x, y, z

  Then
  dn =  |p_in_x| * ε
      + |p_in_y| * ε
      + |p_in_z| * ε
      + |p_out_x| * ε
      + |p_out_y| * ε
      + |p_out_z| * ε

  Compute uncertainty in ψ using propagation of error formula:
  |δψ| = |dn/(sqrt(1-n^2))|
  |δψ| = |-ε/(sqrt(1-n^2))| (|p_in_x| + |p_in_y| + |p_in_z| + |p_out_x| + |p_out_y| + |p_out_z|)

  or more conservative estimate:

  |δψ| = |-ε/(sqrt(1-n^2))| sqrt((1 + |p_in_x|^2 + |p_in_y|^2 + |p_in_z|^2))
  """

  return (180.0/numpy.pi)*(ε/(1+n**2))*(numpy.sum(numpy.abs(p_in)) + numpy.sum(numpy.abs(p_out), axis=1))


def _print_and_write(xform, dfs, dir_table, libs_avail):

  df_info = {
    'values': ['Δθ', 'in degrees'],
    'uncert': ['Δθ-uncert', '(uncert in Δθ in degrees)'],
    'diffs': ['Δθ-diff', 'Δθ - (Δθ geopack_08_dp) in degrees']
  }

  for key in df_info.keys():
    df_str = dfs[xform][key].to_string()
    print('-' * os.get_terminal_size().columns)
    log.info(f"\n{xform} {df_info[key][1]}\n{df_str}")
    file_table = os.path.join(dir_table, f'{xform}_{df_info[key][0]}.txt')
    print(f"Writing above table to {file_table}")
    utilrsw.write(file_table, df_str)

    print("\nStatistics for previous table:")
    print(dfs[xform][key].describe())
    print('-' * os.get_terminal_size().columns)
    print("")


def _config():

  import sys
  from config import CFG

  run_id = 1

  frames = ['GEO', 'MAG', 'GSE', 'GSM']

  if len(sys.argv) > 1:
    run_id = int(sys.argv[1])

  if run_id == 1:
    # Paper run
    axis = 'z'
    delta = {'days': 1}
    to = datetime.datetime(2010, 1, 1, 0, 0, 0)
    tf = datetime.datetime(2015, 1, 1, 0, 0, 0)
    excludes = ['sscweb', 'cxform']

  if run_id == 2:
    # Short run
    axis = 'z'
    delta = {'days': 1}
    to = datetime.datetime(2010, 1, 1, 0, 0, 0)
    tf = datetime.datetime(2010, 1, 3, 0, 0, 0)
    excludes = []

  if run_id == 3:
    # Short high-cadence run
    axis = 'z'
    delta = {'minutes': 10}
    to = datetime.datetime(2010, 12, 21, 0, 0, 0)
    tf = datetime.datetime(2010, 12, 23, 0, 0, 0)
    excludes = ['cxform']

  to_str = datetime.datetime.strftime(to, '%Y%m%d')
  tf_str = datetime.datetime.strftime(tf, '%Y%m%d')
  delta_unit = list(delta.keys())[0]
  delta_str = f'{delta[delta_unit]}{delta_unit}'

  run_str = f'{axis}-delta={delta_str}_{to_str}-{tf_str}'
  dir_table = os.path.join(CFG['cache_dir'], 'angles', 'figures', run_str)
  file_out = os.path.join(CFG['cache_dir'], 'angles', f'{run_str}.pkl')

  print(f'Run id = {run_id}')

  return {
    "axis": axis,
    "delta": delta,
    "to": to,
    "tf": tf,
    "excludes": excludes,
    "frames": frames,
    "dir_table": dir_table,
    "file_out": file_out
  }


cfg = _config()

dfs = {}
for combination in list(itertools.combinations(cfg['frames'], 2)):

  frame_in, frame_out = combination

  hxform_args = {
    'ctype_in': 'car',
    'ctype_out': 'car',
    'frame_in': frame_in,
    'frame_out': frame_out
  }

  libs_avail = libs(frame_in, frame_out, excludes=cfg['excludes'])

  args = (cfg['to'], cfg['tf'], cfg['axis'], cfg['delta'], libs_avail, hxform_args)

  t_dts, Δθ, δψ = angles(*args)

  xform = f"{frame_in}_{frame_out}"

  dfs[xform] = {
    'values': None,
    'diffs': None,
    'uncert': None
  }

  index = pandas.to_datetime(t_dts)
  dfs[xform]['values'] = pandas.DataFrame(Δθ, index=index, columns=libs_avail)
  dfs[xform]['uncert'] = pandas.DataFrame(δψ, index=index, columns=libs_avail)

  dfs[xform]['diffs'] = diffs(dfs, xform, index, libs_avail)

  _print_and_write(xform, dfs, cfg['dir_table'], libs_avail)


print(f"Writing {cfg['file_out']}")
utilrsw.write(cfg['file_out'], dfs)
