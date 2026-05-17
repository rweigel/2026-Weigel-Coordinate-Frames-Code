import os
import numpy
import utilrsw

axis = 'z'
run = f'{axis}-delta=1days_20100101-20150101'

from config import CFG
in_file = os.path.join(CFG['cache_dir'],'angles', f'{run}.pkl')
out_dir = os.path.join(CFG['cache_dir'],'angles', 'figures', run)


def fig_prep():
  from matplotlib import pyplot as plt
  # https://onlinelibrary.wiley.com/page/journal/21699402/homepage/graphics.htm
  plt.gcf().set_size_inches(170/25.4, 228/25.4)
  gs = plt.gcf().add_gridspec(3, hspace=0.07)
  axes = gs.subplots(sharex=True)
  return axes


def plot(df, tranform_str):

  from matplotlib.dates import DateFormatter, YearLocator

  fontsize = 18

  line_map = {
    'geopack_08_dp': ['black', '-'],
    'spacepy': ['blue', '-'],
    'spacepy-irbem': ['blue', '--'],
    'spiceypy1': ['red', '-'],
    'spiceypy2': ['red', '--'],
    'sunpy': ['orange', '-'],
    'pyspedas': ['green', '-'],
    'sscweb': ['purple', '-'],
    'cxform': ['brown', '-'],
    '|max-min|': ['black', '-']
  }

  axes = fig_prep()

  lib = 'geopack_08_dp'
  kwargs = {
    'label': lib,
    'color': line_map[lib][0],
    'linestyle': line_map[lib][1]
  }
  # Plot and adjust y-limits and major ticks for axes[0]
  axes[0].plot(df['values'].index, df['values'][lib], **kwargs)
  axes[0].grid(True)
  axes[0].set_ylabel(f"{tranform_str} [deg]")

  for column in df['diffs'].columns:
    if column == '|max-min|':
      continue

    mabs = utilrsw.mpl.format_exponent(numpy.mean(numpy.abs(df['diffs'][column])), 0)
    label = f"{column} (${mabs}$)"
    kwargs = {
      'label': label,
      'color': line_map[column][0],
      'linestyle': line_map[column][1]
    }
    axes[1].plot(df['diffs'].index, df['diffs'][column], **kwargs)

  axes[1].grid(True)
  axes[1].set_ylabel('$\\Delta$ [deg]')

  # Add zero line to the difference subplot
  axes[1].axhline(0, color='black', linestyle='-', linewidth=1, zorder=0)

  # Force symmetric y-limits for the difference subplot
  yl = axes[1].get_ylim()
  ymax = abs(max(yl, key=abs))
  axes[1].set_ylim(-ymax, ymax)

  axes[1].grid(which='minor', axis='y', linestyle=':', linewidth=0.5)
  kwargs = {
    "ncols": 2,
    "columnspacing": 0.65,
    "handlelength": 1.0,
    "borderaxespad": 0.0,
    "loc": 'upper center'
  }
  axes[1].legend(**kwargs)

  # Plot and adjust y-limits and major ticks for axes[2]
  kwargs = {
    'color': line_map['|max-min|'][0],
    'linestyle': line_map['|max-min|'][1]
  }
  axes[2].plot(df['diffs'].index, df['diffs']['|max-min|'], **kwargs)

  axes[2].grid(True)
  axes[2].set_ylabel('|max-min| [deg]')
  axes[2].set_xlabel('Year')
  axes[2].xaxis.set_major_locator(YearLocator())
  axes[2].xaxis.set_major_formatter(DateFormatter('%Y'))
  yl = axes[2].get_ylim()
  axes[2].set_ylim(0, yl[1])

  utilrsw.mpl.set_fontsize(axes, fontsize=fontsize)
  utilrsw.mpl.adjust_axes(axes)
  utilrsw.mpl.adjust_legend(axes)

  fig = axes[0].get_figure()
  fig.align_ylabels()


utilrsw.mpl.plt_config()
data = utilrsw.read(in_file)

for transform_key in list(data.keys()):
  df = data[transform_key]
  frames = transform_key.split('_')
  frame1 = frames[0]
  frame2 = frames[1]
  axis = axis.upper()
  pair = f"(${axis}_{{{frame1}}}$, ${axis}_{{{frame2}}}$)"
  tranform_str = fr"$\angle$ {pair}"

  print(f"Plotting {transform_key}")
  plot(df, tranform_str)

  utilrsw.mpl.savefig(f'{transform_key}', fdir=out_dir, subdirs=['svg', 'png'], bbox_inches='tight')
