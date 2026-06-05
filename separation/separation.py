import os
import sys

import pandas
import numpy as np

import hxform

from hapiclient import hapi
from hapiclient import hapitime2datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import CFG

R_E = CFG['R_E']

def report(sc, times, positions, title, fname):

  for i in range(4):
    r = np.linalg.norm(positions[i], axis=1)
    r_ave = np.mean(r)
    hxform.xprint(f'  {sc}{i+1} mean distance: {r_ave:.2f} km')
    centroid = np.mean(np.array(positions), axis=0)
    min_dis = np.min(np.linalg.norm(positions[i]-centroid, axis=1))
    hxform.xprint(f'    Minimum distance from centroid: {min_dis:.2f} km')

  separations = []
  labels = []
  for i in range(4):
    for j in range(i+1, 4):
      label = f'{sc}{i+1}-{sc}{j+1}'
      labels.append(label)

      separation = np.linalg.norm(positions[i] - positions[j], axis=1)
      separations.append(separation)

      min_idx = np.argmin(separation)
      min_sep = separation[min_idx]
      min_time = times[i][min_idx]
      if getattr(min_time, 'tzinfo', None) is not None:
        min_time = min_time.replace(tzinfo=None)
      min_time_str = min_time.isoformat()
      hxform.xprint(f"  {label}")
      angle = np.degrees(min_sep/r_ave)
      hxform.xprint(f'    Minimum separation distance and angle: {min_sep:.2f} km and {angle:.2e} deg at {min_time_str}')

  plot(sc, times, separations, labels, title, fname)


def plot(sc, times, separations, labels, title, fname):
  from matplotlib import pyplot as plt
  from datetick import datetick

  for time, separation, label in zip(times, separations, labels):
    plt.plot(time, separation, label=label)

  plt.title(title)
  plt.ylabel('Separation Distance (km)')
  plt.grid()
  plt.legend()
  datetick('x')

  writefig(sc, fname)
  plt.close()


def writefig(sc, fname):
  from matplotlib import pyplot as plt
  import os

  for fmt in ['pdf', 'svg', 'png']:
    subdir = fmt
    if fmt == 'pdf':
      subdir = ''
    fname_full = os.path.join(CFG['cache_dir'], 'separation', sc.lower(), 'figures', subdir, fname)

    if not os.path.exists(os.path.dirname(fname_full)):
      os.makedirs(os.path.dirname(fname_full))

    hxform.xprint(f'  Writing: {os.path.relpath(fname_full)}.{fmt}')
    if fmt == 'png':
      plt.savefig(f'{fname_full}', dpi=300)
    else:
      plt.savefig(f'{fname_full}.{fmt}')


opts  = {
  'logging': CFG['hapi_logging'],
  'usecache': True,
  'cachedir': os.path.join(CFG['cache_dir'], 'positions', 'hapi')
}

def _cdaweb(sc, start, stop, frame):

  server = 'https://cdaweb.gsfc.nasa.gov/hapi'

  if sc == 'MMS':
    dataset_suffix = 'MEC_SRVY_L2_EPHT89D'
    parameters_suffix = f'mec_r_{frame.lower()}'

    # EPD data has time stamps that differ between s/c
    #dataset = f'{sc.upper()}_EPD-EIS_SRVY_L2_ELECTRONENERGY'
    #parameters = f'{sc}_epd_eis_srvy_l2_electronenergy_position_{frame_cdaweb.lower()}'

  if sc == 'Cluster':
    dataset_suffix = 'CP_FGM_5VPS'
    parameters_suffix = f'cp_r_{frame.lower()}'

  result = {}

  result['title'] = f'CDAWeb {frame} ({sc} {dataset_suffix})'
  result['fname'] = f'separation_CDAWeb_{frame}_{dataset_suffix}'

  hxform.xprint('')
  hxform.xprint(result['title'])


  times = []
  positions = []
  dfs = []

  for s in range(1, 5):
    if sc == 'MMS':
      dataset = f'{sc}{s}_{dataset_suffix}'
      parameters = f'{sc.lower()}{s}_{parameters_suffix}'
    if sc == 'Cluster':
      dataset = f'C{s}_{dataset_suffix}'
      parameters = f'sc_pos_xyz_{frame.lower()}__{dataset}'

    print(f"  Getting CDAWeb {dataset} {parameters}...")

    data, meta = hapi(server, dataset, parameters, start, stop, **opts)
    times.append(hapitime2datetime(data['Time']))
    positions.append(data[parameters])
    dfs.append(pandas.DataFrame(data[parameters], columns=['x', 'y', 'z'], index=times[-1]))

  # Keep only timestamps present in all four spacecraft so arrays are aligned.
  time_to_idx = [{t: j for j, t in enumerate(t_arr)} for t_arr in times]
  common_times = [t for t in times[0] if all(t in m for m in time_to_idx[1:])]

  if len(common_times) == 0:
    raise ValueError("No common timestamps across spacecraft for CDAWeb data.")
  else:
    hxform.xprint(f"  Found {len(common_times)}/{len(times[0])} common timestamps across spacecraft for CDAWeb data.")

  for i in range(4):
    idx = [time_to_idx[i][t] for t in common_times]
    p_arr = np.asarray(positions[i])
    times[i] = np.asarray(common_times)
    positions[i] = p_arr[idx]

  result['times'] = times
  result['positions'] = positions
  result['dfs'] = dfs

  return result


def _sscweb(sc, start, stop, frame):
  server = 'http://hapi-server.org/servers/SSCWeb/hapi'

  result = {}
  result['title'] = f'SSCWeb {frame}'
  result['fname'] = f'separation_SSCWeb_{frame}'

  hxform.xprint('')
  hxform.xprint(result['title'])

  times = []
  positions = []
  dfs = []

  if sc == 'Cluster':
    sc = sc.lower()

  for s in range(1, 5):
    dataset    = f'{sc}{s}'
    parameters = f'X_{frame},Y_{frame},Z_{frame}'

    print(f"  Getting SSCWeb {dataset} {parameters}...")
    data, meta = hapi(server, dataset, parameters, start, stop, **opts)

    xyz = np.column_stack((data[f'X_{frame}'], data[f'Y_{frame}'], data[f'Z_{frame}']))

    times.append(hapitime2datetime(data['Time']))

    positions.append(xyz*R_E)
    dfs.append(pandas.DataFrame(xyz*R_E, columns=['x', 'y', 'z'], index=times[-1]))

  result['times'] = times
  result['positions'] = positions
  result['dfs'] = dfs

  return result

sc = 'Cluster'
#sc = 'MMS'

if sc == 'MMS':
  start = '2016-09-14'
  stop  = '2016-09-16'

  # "Its four spacecraft are flying only four-and-a-half miles apart" => 7.24 km
  angle = 7.24/(8.79*CFG['R_E'])
  hxform.xprint(f"4.5 miles ({7.24} km) separation => {np.degrees(angle):.2e} deg at 8.79 R_E")

  """
  https://www.nasa.gov/missions/mms/nasas-mms-achieves-closest-ever-flying-formation/
  On Sept. 15, 2016, NASA’s Magnetospheric Multiscale, or MMS, mission achieved
  a new record: Its four spacecraft are flying only four-and-a-half miles apart,
  the closest separation ever of any multi-spacecraft formation.

  In the following, we compute the separation distances and angles from s/c
  pairs using position data from CDAWeb and SSCWeb. "closest separation" seems
  to mean the minimum distance of any s/c from centroid of the four s/c.

  See also https://mmsvis.gsfc.nasa.gov/ for plots
  """

if sc == 'Cluster':
  """
  https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2021JA029474
  "With careful maneuvering and monitoring the team managed to achieve a
  separation of 4 km between Cluster 3 and Cluster 4 on September 19, 2013.
  This firstachievement paved the way to future tight formations to study the
  Earth's bow shock at even smaller scales: with 3 km separation reached in
  January 2016 and 2.5 km achieved in December 2018. This is the smallest
  distance ever achieved between two magnetospheric spacecraft without active
  ranging or GPS receivers.
  """

  #start = '2013-09-19'
  #stop  = '2013-09-20'
  start = '2016-01-01'
  stop  = '2016-02-01'
  start = '2018-12-01'
  stop  = '2019-01-01'

hxform.xprint(f"Start: {start}, Stop: {stop}")

#frame = 'GSM'
#frame_cdaweb = 'ECI'
frame_cdaweb = 'GSE'
if frame_cdaweb == 'ECI':
  frame_sscweb = 'J2K'
else:
  frame_sscweb = frame_cdaweb


cdaweb = _cdaweb(sc, start, stop, frame_cdaweb)
report(sc, cdaweb['times'], cdaweb['positions'], cdaweb['title'], cdaweb['fname'])

sscweb = _sscweb(sc, start, stop, frame_sscweb)
report(sc, sscweb['times'], sscweb['positions'], sscweb['title'], sscweb['fname'])


hxform.xprint('')
hxform.xprint('First few timestamps and positions for each source:')

for i in range(4):
  hxform.xprint('')
  hxform.xprint(40*'-')
  hxform.xprint(f'{sc}{i+1}')
  hxform.xprint(40*'-')

  # Show first few lines of each dataframe
  hxform.xprint(f"CDAWeb {frame_cdaweb}:")
  hxform.xprint(cdaweb['dfs'][i].head())

  hxform.xprint(f"\nSSCWeb {frame_sscweb}:")
  hxform.xprint(sscweb['dfs'][i].head())


# Move separation.log to fname.log
fname_log = os.path.join(CFG['cache_dir'], 'separation', sc.lower(), 'separation.log')
os.rename('separation.log', fname_log)


