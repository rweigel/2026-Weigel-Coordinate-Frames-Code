import os
import sys
import pandas
import numpy as np

import hxform

from hapiclient import hapi
from hapiclient import hapitime2datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import CFG

# https://www.nasa.gov/missions/mms/nasas-mms-achieves-closest-ever-flying-formation/
# On Sept. 15, 2016, NASA’s Magnetospheric Multiscale, or MMS, mission achieved
# a new record: Its four spacecraft are flying only four-and-a-half miles apart,
# the closest separation ever of any multi-spacecraft formation.

# In the following, we compute the separation distances and angles from s/c
# pairs using position data from CDAWeb and SSCWeb. "closest separation" seems
# to mean the minimum distance of any s/c from centroid of the four s/c.

# See also https://mmsvis.gsfc.nasa.gov/ for plots

# "Its four spacecraft are flying only four-and-a-half miles apart" => 7.24 km
angle = 7.24/(8.79*CFG['R_E'])

hxform.xprint(f"4.5 miles ({7.24} km) separation => {np.degrees(angle):.2e} deg at 8.79 R_E")

R_E = CFG['R_E']

def report(times, positions, title, fname):

  for i in range(4):
    r = np.linalg.norm(positions[i], axis=1)
    r_ave = np.mean(r)
    hxform.xprint(f'  MMS{i+1} mean distance: {r_ave/R_E:.2f} R_E')
    centroid = np.mean(np.array(positions), axis=0)
    min_dis = np.min(np.linalg.norm(positions[i]-centroid, axis=1))
    hxform.xprint(f'    Minimum distance from centroid: {min_dis:.2f} km')

  separations = []
  labels = []
  for i in range(4):
    for j in range(i+1, 4):
      label = f'MMS{i+1}-MMS{j+1}'
      labels.append(label)

      separation = np.linalg.norm(positions[i] - positions[j], axis=1)
      separations.append(separation)

      min_sep = np.min(separation)
      hxform.xprint(f"  {label}")
      angle = np.degrees(min_sep/r_ave)
      hxform.xprint(f'    Minimum separation distance and angle: {min_sep:.2f} km and {angle:.2e} deg')

  plot(times, separations, labels, title, fname)

def plot(times, separations, labels, title, fname):
  from matplotlib import pyplot as plt
  from hapiplot.plot.datetick import datetick

  for time, separation, label in zip(times, separations, labels):
    plt.plot(time, separation, label=label)

  plt.title(title)
  plt.ylabel('Separation Distance (km)')
  plt.grid()
  plt.legend()
  datetick(dir='x')

  writefig(fname)
  plt.close()

def writefig(fname):
  from matplotlib import pyplot as plt
  import os

  for fmt in ['pdf', 'svg', 'png']:
    subdir = fmt
    if fmt == 'pdf':
      subdir = ''
    fname_full = os.path.join(CFG['cache_dir'], 'mms_separation', 'figures', subdir, fname)

    if not os.path.exists(os.path.dirname(fname_full)):
      os.makedirs(os.path.dirname(fname_full))

    hxform.xprint(f'  Writing: {fname_full}.{fmt}')
    if fmt == 'png':
      plt.savefig(f'{fname_full}', dpi=300)
    else:
      plt.savefig(f'{fname_full}.{fmt}')


start = '2016-09-14'
stop  = '2016-09-16'
opts  = {
  'logging': CFG['hapi_logging'],
  'usecache': True,
  'cachedir': os.path.join(CFG['cache_dir'], 'positions', 'hapi')
}

print("Start: {start}, Stop: {stop}".format(start=start, stop=stop))

#frame = 'GSM'
frame_cdaweb = 'GSE'
#frame_cdaweb = 'ECI'

server = 'https://cdaweb.gsfc.nasa.gov/hapi'
dataset_suffix = 'MEC_SRVY_L2_EPHT89D'
title = f'CDAWeb {frame_cdaweb} (MMSi_{dataset_suffix})'
fname = 'mms_separation_' + title.replace(" ", "_").replace("(", "").replace(")", "")

hxform.xprint('')
hxform.xprint(title)

times = {'cdaweb': [], 'sscweb': []}
positions = {'cdaweb': [], 'sscweb': []}
dfs = {'cdaweb': [], 'sscweb': []}

for s in range(1, 5):
  sc   = f'mms{s}'
  dataset = f'{sc.upper()}_MEC_SRVY_L2_EPHT89D'
  parameters = f'{sc}_mec_r_{frame_cdaweb.lower()}'
  # EPD data has time stamps that differ between s/c
  #dataset = f'{sc.upper()}_EPD-EIS_SRVY_L2_ELECTRONENERGY'
  #parameters = f'{sc}_epd_eis_srvy_l2_electronenergy_position_{frame_cdaweb.lower()}'
  print(f"  Getting CDAWeb {dataset} {parameters}...")
  data, meta = hapi(server, dataset, parameters, start, stop, **opts)
  times['cdaweb'].append(hapitime2datetime(data['Time']))
  positions['cdaweb'].append(data[parameters])
  dfs['cdaweb'].append(pandas.DataFrame(data[parameters], columns=['x', 'y', 'z'], index=times['cdaweb'][-1]))


report(times['cdaweb'], positions['cdaweb'], title, fname)

server = 'http://hapi-server.org/servers/SSCWeb/hapi'
if frame_cdaweb == 'ECI':
  frame_sscweb = 'J2K'
else:
  frame_sscweb = frame_cdaweb

title = f'SSCWeb {frame_sscweb}'

hxform.xprint('')
hxform.xprint(title)

for s in range(1, 5):
  dataset    = f'mms{s}'
  parameters = f'X_{frame_sscweb},Y_{frame_sscweb},Z_{frame_sscweb}'

  print(f"  Getting SSCWeb {dataset} {parameters}...")
  data, meta = hapi(server, dataset, parameters, start, stop, **opts)

  xyz = np.column_stack((data[f'X_{frame_sscweb}'], data[f'Y_{frame_sscweb}'], data[f'Z_{frame_sscweb}']))

  times['sscweb'].append(hapitime2datetime(data['Time']))
  positions['sscweb'].append(xyz*R_E)
  dfs['sscweb'].append(pandas.DataFrame(xyz*R_E, columns=['x', 'y', 'z'], index=times['sscweb'][-1]))

fname = 'mms_separation_' + title.replace(" ", "_").replace("(", "").replace(")", "")

report(times['sscweb'], positions['sscweb'], title, fname)

hxform.xprint('')
hxform.xprint('First few timestamps and positions for each source:')

for i in range(4):
  hxform.xprint('')
  hxform.xprint(40*'-')
  hxform.xprint(f'MMS{i+1}')
  hxform.xprint(40*'-')

  # Show first few lines of each dataframe
  hxform.xprint(f"CDAWeb {frame_cdaweb} {dataset_suffix}:")
  hxform.xprint(dfs['cdaweb'][i].head())

  hxform.xprint(f"\nSSCWeb {frame_sscweb}:")
  hxform.xprint(dfs['sscweb'][i].head())
