
import utilrsw

from config import CFG


def _compute_diffs(info1, info2, opts):
  import numpy as np

  from collections import Counter

  print(f"    Computing differences between {info1['name']} and {info2['name']}")

  info1['time'] = np.array([t.replace(tzinfo=None) for t in info1['time']])
  info2['time'] = np.array([t.replace(tzinfo=None) for t in info2['time']])

  def _interp_dts(times1, times2):
    # For each time in time1, compute minimum absolute difference to any time in time2
    dt_mins = np.empty(len(times1), dtype=float)
    for i, time1 in enumerate(times1):
      dt_mins[i] = np.min(np.abs(time1-times2)).total_seconds()
    return dt_mins

  dt1_mins = _interp_dts(info1['time'], info2['time'])
  print(f"      {info1['name']} max dt to {info2['name']}: {np.max(dt1_mins):.3f} s")
  dt2_mins = _interp_dts(info2['time'], info1['time'])
  print(f"      {info2['name']} max dt to {info1['name']}: {np.max(dt2_mins):.3f} s")

  # Compute time differences in seconds between consecutive samples in each source
  dt1 = np.diff([t.timestamp() for t in info1['time']])
  dt2 = np.diff([t.timestamp() for t in info2['time']])

  # Round each dt to 3 decimal places
  dt1 = np.round(dt1, 3)
  dt2 = np.round(dt2, 3)

  # Compute histogram (as Counter) for each
  hist1 = Counter(dt1)
  hist2 = Counter(dt2)

  # Sort histogram by frequency (descending)
  hist1 = dict(sorted(hist1.items(), key=lambda item: item[1], reverse=True))
  hist2 = dict(sorted(hist2.items(), key=lambda item: item[1], reverse=True))

  # Keep first 10 most common time differences
  hist1 = dict(list(hist1.items())[:10])
  hist2 = dict(list(hist2.items())[:10])

  # Convert keys to float so numpy.float64(value) not printed in output.
  hist1 = {float(k): v for k, v in hist1.items()}
  hist2 = {float(k): v for k, v in hist2.items()}
  print(f"      {info1['name']} dt histogram (first 10) (s): {dict(hist1)}")
  print(f"      {info2['name']} dt histogram (first 10) (s): {dict(hist2)}")

  if opts['interp_times']:
    # Only keep times where min dt to nearest time in other source is less
    # than 10 seconds
    print(f"      Filtering timestamps with dt to nearest time in other source >= {opts['interp_times_tolerance']} seconds.")
    mask1 = dt1_mins < opts['interp_times_tolerance']
    time1 = info1['time'][mask1]
    xyz1  = info1['xyz'][mask1]
    print(f"      {info1['name']} has {len(time1)} timestamps after filtering (originally {len(info1['time'])})")
    mask2 = dt2_mins < opts['interp_times_tolerance']
    time2 = info2['time'][mask2]
    xyz2 = info2['xyz'][mask2]
    print(f"      {info2['name']} has {len(time2)} timestamps after filtering (originally {len(info2['time'])})")

    print("      Computing union of timestamps.")
    # Union of timestamps from CDAWeb and SSCWeb
    union_times = set(time1).union(set(time2))
    union_times = sorted(union_times)

    if len(union_times) == 0:
      print("    No timestamps to compare after filtering. Returning empty arrays.")
      return [], [], [], []

    print(f"      Interpolating to {len(union_times)} timestamps.")
    info1['xyz_i'] = _interp(union_times, time1, xyz1)
    info1['time_i'] = union_times
    info2['xyz_i'] = _interp(union_times, time2, xyz2)
    info2['time_i'] = union_times

    time_key = 'time_i'
    xyz_key = 'xyz_i'
  else:
    # Common timestamps
    print("      Computing common timestamps.")
    common_times = set(info2['time']).intersection(set(info1['time']))
    common_times = sorted(common_times)
    common_times = sorted(common_times)
    print(f"      {len(common_times)} common timestamps.")

    # Find the indices of common times
    idx1 = np.nonzero(np.isin(info1['time'], common_times))[0]
    idx2 = np.nonzero(np.isin(info2['time'], common_times))[0]

    info1['xyz_c'] = np.array(info1['xyz'][idx1])
    info1['time_c'] = common_times
    info2['xyz_c'] = np.array(info2['xyz'][idx2])
    info2['time_c'] = common_times
    time_key = 'time_c'
    xyz_key = 'xyz_c'

    if not opts['interp_times'] and opts['print_common_times']:
      for i in range(len(common_times)):
        t_str = common_times[i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
        print(f"      {t_str} | {info1['name']}: {info1[xyz_key][i]} | {info2['name']}: {info2[xyz_key][i]}")

  nt = len(info1[time_key]) # Number of timestamps
  Δθ = np.full(nt, np.nan)
  Δr = np.full(nt, np.nan)
  t = np.empty(nt, dtype='datetime64[ns]')
  t[:] = np.datetime64('NaT')
  # Average r
  r_ave = 0.5*(np.linalg.norm(info1[xyz_key], axis=1) + np.linalg.norm(info2[xyz_key], axis=1))

  for i in range(nt):
    t[i] = np.datetime64(info1[time_key][i])
    Δr[i] = np.linalg.norm(info1[xyz_key][i] - info2[xyz_key][i])
    # d = denominator for angle calculation
    # Δθ[i] = arccos [ (a·b)/(|a|*|b|) ] = arccos (n/d)
    n = np.dot(info1[xyz_key][i], info2[xyz_key][i])
    d = np.linalg.norm(info1[xyz_key][i])*np.linalg.norm(info2[xyz_key][i])
    if np.abs(n) > np.abs(d):
      if np.abs(n) - np.abs(d) > 1e-10:
        raise ValueError(f"|n| - |d| > 1e-10. n = {n:.16f} d = {d:.16f}, n-d = {n-d:.16f}")
      if opts['cos_warnings']:
        # If n == d, then the angle is 0 degrees, but numpy gives
        # RuntimeWarning: invalid value encountered in arccos. Did rounding
        # cause n > d?
        # TODO: Consider https://github.com/sunpy/sunpy/pull/7530#issuecomment-2020890282
        wmsg =  f"      Warning: {t[i]}:\n"
        wmsg +=  "        |n| > |d| in Δθ[i] = arccos [ (a·b)/(|a|*|b|) ] = arccos (n/d)\n"
        wmsg += f"        a = {info1[xyz_key][i]}\n"
        wmsg += f"        b = {info2[xyz_key][i]}\n"
        wmsg += f"        n = {n:.16f}\n"
        wmsg += f"        d = {d:.16f}\n"
        wmsg += f"        n/d = {n/d}\n"
        wmsg +=  "        Using Δθ = 0\n"
        print(wmsg)
      Δθ[i] = 0
    else:
      Δθ[i] = (180/np.pi)*np.arccos(n/d)

    if opts['interp_times'] and opts['print_interp_vals']:
      t_str = info1[time_key][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
      print(f"      t = {t_str} | Δr = {Δr[i]:.5f} [R_E] | Δ∠: {Δθ[i]:.5f}°")

  return t, r_ave, Δr, Δθ


def _interp(timei, time, xyz):
  from scipy.interpolate import interp1d

  # Interpolate CDAWeb data onto the union of timestamps
  interp_func = interp1d(
    [t.timestamp() for t in time],
    xyz,
    axis=0,
    bounds_error=False,
    fill_value="extrapolate"
  )

  return interp_func([t.timestamp() for t in timei])


def _plot_xyz(ax, info1, info2, t, r_ave):

  print(f"  Plotting xyz for {info1['name']} and {info2['name']}")

  colors = ['r', 'g', 'b']
  component = ['X', 'Y', 'Z']

  for c in range(3):
    label1 = f"{info1['name']}/{component[c]}"
    label2 = f"{info2['name']}/{component[c]}"
    ax.plot(info1['time'], info1['xyz'][:,c],
            label=label1, lw=2, linestyle='-', color=colors[c])
    ax.plot(info2['time'], info2['xyz'][:,c],
            label=label2, lw=3, linestyle='--', color=colors[c])

  #label = '$\\overline{r}$'
  r_ave_line, = ax.plot(t, r_ave, lw=2, linestyle='-', color='k')

  #_adjust_y_range(ax, gap_fraction=1)
  ax.set_ylabel('$R_E$', rotation=0)
  ax.grid()

  lkw = {**CFG['plot_opts']['legend_kwargs'], 'ncols': 3}
  ax.legend(**lkw)

  # Short black line + textbox label for r_ave in the lower left
  ax.plot([0.02, 0.06], [0.08, 0.08], color='k', lw=2,
          transform=ax.transAxes, clip_on=False, solid_capstyle='butt')
  ax.annotate('$\\overline{r}$',
              xy=(0.07, 0.08), xycoords='axes fraction',
              va='center', ha='left',
              fontsize=CFG['plot_opts']['fontsize'],
              bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none'))

  ax.set_xticklabels([])


def _plot_diffs(ax, t, r_ave, Δr, Δθ, R_E):
  import numpy as np
  print("  Plotting diffs")

  lw = 2  # line width

  Δr_rel = Δr/r_ave
  Δr_rel_max = np.nanmax(Δr_rel)
  Δr_rel_max_str = f" (max= 1/{1/Δr_rel_max:.0f})"

  Δr_max = np.nanmax(Δr)
  Δr_max_str = f'($|Δ\\mathbf{{r}}|_{{\\rm{{max}}}} = {Δr_max*R_E:.1f}$ [km])'

  print(f"    Δr_max = {Δr_max:.5f} [R_E]")
  print(f"    Δr_max = {Δr_max*R_E:.1f} [km]")

  ax.plot(t, Δr, 'r-', lw=lw, label=f'$|Δ\\mathbf{{r}}|/R_E$ {Δr_max_str}')
  ax.plot(t, Δr_rel, 'b-', lw=lw, label=f'$|Δ\\mathbf{{r}}|/\\overline{{r}}$ {Δr_rel_max_str}')
  ax.plot(t, Δθ, 'g-', lw=lw, label='$Δθ$ [deg]')

  #_adjust_y_range(ax, bottom=0, gap_fraction=1)
  lkw = {**CFG['plot_opts']['legend_kwargs'], 'ncols': 2, 'handlelength': 1.5}
  ax.legend(**lkw)
  ax.grid()


def _figprep():
  from matplotlib import pyplot as plt

  plt.gcf().set_size_inches(CFG['plot_opts']['figsize_inches'])
  gs = plt.gcf().add_gridspec(2)
  axes = gs.subplots(sharex=True)

  for ax in axes:
    ax.grid(which='minor', linestyle=':', linewidth=0.5, color=[0.75]*3)
    ax.minorticks_on()

  return axes


def _savefigs(fname):
  import os
  from matplotlib import pyplot as plt

  for fmt in ['svg', 'png', 'pdf']:
    bbox_inches = None
    #kwargs = {'bbox_inches': 'tight'}
    kwargs = {'bbox_inches': bbox_inches}
    if fmt == 'png':
      kwargs['dpi'] = 300

    base = f'{CFG["cache_dir"]}/positions/figures'
    if fmt == 'pdf':
      fname_full = f'{base}/{fname}.{fmt}'
    else:
      fname_full = f'{base}/{fmt}/{fname}.{fmt}'
    os.makedirs(os.path.dirname(fname_full), exist_ok=True)
    print(f"  Writing {os.path.relpath(fname_full)}")
    plt.savefig(fname_full, bbox_inches=bbox_inches)
  plt.close()


def _plot(satellite, info1, info2, opts):

  from datetick import datetick

  axes = _figprep()

  title = f"{satellite} {info1['name']}/{info1['frame']} "
  title += f"vs. {info2['name']}/{info2['frame']}"
  #axes[0].set_title(title)

  fname = f"{satellite}_{info1['name']}-{info1['frame']}"
  fname += f"_vs_{info2['name']}-{info2['frame']}"

  t, r_ave, Δr, Δθ = _compute_diffs(info1, info2, opts)

  if len(t) == 0:
    print("  No data to plot after filtering. Skipping plot.")
    return

  _plot_xyz(axes[0], info1, info2, t, r_ave)
  _plot_diffs(axes[1], t, r_ave, Δr, Δθ, opts['R_E'])

  utilrsw.mpl.set_fontsize(axes, CFG['plot_opts']['fontsize'])
  utilrsw.mpl.adjust_axes(axes)
  utilrsw.mpl.adjust_legend(axes, debug=True)

  datetick('x', adjust_first_xlabel=True, adjust_last_xlabel=True)

  _savefigs(fname)


def plot_combos(satellite, sources, opts):
  utilrsw.mpl.plt_config()

  for i, source1 in enumerate(sources):
    for j, source2 in enumerate(sources):
      if i >= j or source1 is None or source2 is None:
        continue
      s1 = f"{source1['name']}/{source1['frame']}"
      s2 = f"{source2['name']}/{source2['frame']}"
      print(f"  Comparing {s1} with {s2}")
      _plot(satellite, source1, source2, opts)


if __name__ == "__main__":
  from config import cli
  args = cli()

  if args.satellites == CFG['known_satellites']:
    base = f"{CFG['cache_dir']}/positions/positions"
  else:
    base = f"{CFG['cache_dir']}/positions/positions-{'_'.join(args.satellites)}"

  print(f"Reading position data from {base}.pkl")
  try:
    results = utilrsw.read(f'{base}.pkl')
  except:
    if args.satellites == CFG['known_satellites']:
      msg = f"No cached data found at {base}.pkl. Run python positions.py to generate the data before plotting."
      print(msg)
    else:
      msg = f"No cached data found at {base}.pkl. Run python positions.py --satellites "
      msg += f"{','.join(args.satellites)} to generate the data before plotting."
      print(msg)
    exit(1)

  for satellite in args.satellites:
    for frame in results[satellite]:
      sources = results[satellite][frame]
      print(f"Plotting position comparisons for {satellite}/{frame}")
      plot_combos(satellite, sources, CFG['plot_opts'])
