import os

import numpy as np

import utilrsw

from hapiclient import hapi
from hapiclient import hapitime2datetime

import warnings
warnings.filterwarnings("ignore", message="The argument 'infer_datetime_format'")
warnings.filterwarnings("ignore", message="Could not infer format")

# Set to ERROR to suppress INFO message for requests to Horizons server
import sunpy
sunpy.log.setLevel('ERROR')

from config import CFG, cli

def _infos(satellite):
  """
  Returns dict with keys of data sources
    {
      'sscweb': [info1, info2, ...],
      'cdaweb': [info1, info2, ...],
      'amda': [info1, info2, ...],
      'jpl': [info1, info2, ...]
    }
  where each info is a dict with info needed to request data from the respective
  source. Each info subscript corresponds to a different frame. If a data source
  does not have data for a given frame, then the corresponding info is None.

  Comparisons are made between data obtained using info1 for all data sources, then
  info2 for all data sources, etc.

  """

  if satellite not in CFG['known_satellites']:
    raise ValueError(f"Unknown satellite: {satellite}. Must be one of {CFG['known_satellites']}")

  infos = {'cdaweb': [], 'sscweb': [], 'amda': [], 'jpl': []}

  if satellite.startswith('Cluster'):

    sc = f'cluster{satellite[-1]}'
    dataset = f'C{satellite[-1]}_CP_FGM_5VPS'
    start = '2010-09-01T09:09:00.100Z'
    stop  = '2010-09-02T00:00:00.000Z'

    for frame in ['GSE', 'GSM']:

      if frame == 'GSE':
        infos['cdaweb'].append({
          'name': 'CDAWeb',
          'dataset': dataset,
          'parameter': f'sc_pos_xyz_{frame.lower()}__{dataset}',
          'start': start,
          'stop':  stop,
          'frame': frame
        })
      else:
        # No GSM data in CDAWeb
        infos['cdaweb'].append(None)

      infos['sscweb'].append({
        'name': 'SSCWeb',
        'dataset': sc,
        'start': start,
        'stop':  stop,
        'frame': frame
      })

      infos['jpl'].append(None)

      infos['amda'].append({
        'name': 'AMDA',
        'dataset': f'clust{satellite[-1]}-orb-all',
        'parameter': f'c{satellite[-1]}_xyz_{frame.lower()}',
        'start': start,
        'stop':  stop,
        'frame': frame
      })

  if satellite == 'DSCOVR':
    start = '2021-11-25T00:00:00Z'
    stop  = '2021-12-05T00:00:00Z'

    for frame in ['GSE', 'GCI']: # GCI must be last
      infos["cdaweb"].append({
        'name': 'CDAWeb',
        'dataset': 'DSCOVR_ORBIT_PRE',
        'parameter': f'{frame}_POS',
        'start': start,
        'stop':  stop,
        'frame': frame
      })

      infos["sscweb"].append({
          'name': 'SSCWeb',
          'dataset': 'dscovr',
          'start': start,
          'stop':  stop,
          'frame': frame
        })

      infos["jpl"].append({
        'name': 'JPL',
        'spacecraft': 'dscovr',
        'start': start,
        'stop':  stop,
        'frame': frame if frame != 'GCI' else 'GEI'
      })

      infos["amda"].append(None)

    infos["cdaweb"].append(infos["cdaweb"][-1].copy())
    infos["sscweb"].append(infos["sscweb"][-1].copy())
    infos["sscweb"][-1]['frame'] = 'J2K'
    infos["jpl"].append(infos["jpl"][-1].copy())
    infos['amda'].append(None)

  if satellite.startswith('THEMIS'):

    sc = f'themis{satellite[-1].lower()}'
    start = '2016-09-01T00:00:00Z'
    stop  = '2016-09-02T00:00:00Z'
    for frame in ['GSE', 'GSM', 'GEI']: # GEI must be last
      if frame == 'GEI':
        # Metadata indicates COORDINATE_SYSTEM = GEI for this parameter
        parameter = f'th{sc[-1]}_pos'
      else:
        parameter = f'th{sc[-1]}_pos_{frame.lower()}'

      infos['cdaweb'].append({
        'name': 'CDAWeb',
        'dataset': f'TH{sc[-1].upper()}_L1_STATE@0',
        'parameter': parameter,
        'start': start,
        'stop':  stop,
        'frame': frame
      })

      infos['sscweb'].append({
        'name': 'SSCWeb',
        'dataset': sc,
        'start': start,
        'stop':  stop,
        'frame': frame
      })

      if satellite[-1].lower() in ['b', 'c']:
        infos['jpl'].append({
          'name': 'JPL',
          'spacecraft': sc,
          'start': start,
          'stop':  stop,
          'frame': frame
        })
      else:
        infos['jpl'].append(None)

      infos['amda'].append(None)

    infos['cdaweb'].append(infos['cdaweb'][-1].copy())
    infos['sscweb'].append(infos['sscweb'][-1].copy())
    infos['sscweb'][-1]['frame'] = 'J2K'
    infos['jpl'].append(None)
    infos['amda'].append(None)

  if satellite.startswith('MMS'):

    # TODO:
    # MMS1_MEC_SRVY_L2_EPHT89D
    # has GEI, GSM, GEO, SM
    sc   = f'mms{satellite[-1]}'
    start = '2016-09-01T00:00:00Z'
    stop  = '2016-09-02T00:00:00Z'

    for frame in ['GEI', 'GSE', 'GSM']:

      if frame in ['GSE', 'GSM']:
        infos['cdaweb'].append({
          'name': 'CDAWeb',
          'dataset': f'{sc.upper()}_EPD-EIS_SRVY_L2_ELECTRONENERGY',
          'parameter': f'{sc}_epd_eis_srvy_l2_electronenergy_position_{frame.lower()}',
          'start': start,
          'stop':  stop,
          'frame': frame
        })
      else:
        # No J2K or GEI data in CDAWeb dataset
        infos['cdaweb'].append(None)

      infos['sscweb'].append({
        'name': 'SSCWeb',
        'dataset': sc,
        'start': start,
        'stop':  stop,
        'frame': frame
      })

      infos['jpl'].append({
        'name': 'JPL',
        'spacecraft': sc,
        'start': start,
        'stop':  stop,
        'frame': frame
      })

      infos['amda'].append(None)

  if satellite == 'Geotail':
    """
    https://cdaweb.gsfc.nasa.gov/misc/NotesG.html#GE_OR_DEF
    CDAWeb dataset has GCI, GSE, GSM, HEC
    """

    start      = '2021-11-25T00:00:00Z'
    stop       = '2021-12-05T00:00:00Z'

    for frame in ['GSE', 'GSM', 'GCI']: # GCI must be last
      infos['cdaweb'].append({
        'name': 'CDAWeb',
        'dataset': 'GE_OR_DEF',
        'parameter': f'{frame}_POS',
        'start': start,
        'stop':  stop,
        'frame': frame
      })

      infos['sscweb'].append({
        'name': 'SSCWeb',
        'dataset': 'geotail',
        'start': start,
        'stop':  stop,
        'frame': frame if frame != 'GCI' else 'GEI'
      })

      infos['jpl'].append(None)

      infos['amda'].append(None)

    # Also compare GCI from CDAWeb with J2K from SSCWeb (which is why GCI must
    # be last in the loop above).
    infos['cdaweb'].append(infos['cdaweb'][-1].copy())
    infos['sscweb'].append(infos['sscweb'][-1].copy())
    infos['sscweb'][-1]['frame'] = 'J2K'
    infos['jpl'].append(None)
    infos['amda'].append(None)

  return infos


def _jpl(info, time, logging=False):
  import pickle

  from sunpy.coordinates import get_horizons_coord

  if info is None:
    return None

  if info['frame'] not in ['GSE', 'GEI', 'GSM']:
    raise ValueError(f"Unknown frame: {info['frame']}. Must be one of 'GSE', 'GEI', or 'GSM'.")

  # To find ids, see https://ssd.jpl.nasa.gov/horizons/app.html#/
  # and edit the "Target Body" field to find the id. Does not seem possible
  # to query for list of Target Body ids.
  # Note: No Geotail or themis{a,d}
  known_spacecraft = {
    'ace': -92,
    'dscovr': -78,
    'mms1': -140482,
    'mms2': -140483,
    'mms3': -140484,
    'mms4': -140485,
    'themisb': -192,
    'themisc': -193
  }

  if info['spacecraft'] not in known_spacecraft:
    raise ValueError(f"Unknown target body: {info['target_body']}. Must be one of {list(known_spacecraft.keys())}")

  target_body = known_spacecraft[info['spacecraft']]
  to = info['start']
  tf = info['stop']
  cache_file = f"{CFG['cache_dir']}/positions/jpl/jpl-{info['spacecraft']}-{info['frame']}-{to}-{tf}.pkl"

  if os.path.exists(cache_file):
    print(f"Loading JPL data from cache file: {cache_file}")
    with open(cache_file, 'rb') as f:
      return pickle.load(f)

  #solar_system_positions.set('de432s')
  #solar_system_positions.set('de440s')
  dt = time[1] - time[0]
  # Using, e.g., 60s gives error: https://github.com/sunpy/sunpy/issues/8188
  step_min = int(dt.total_seconds()/60)
  t = {
    "start": time[0],
    "stop": time[-1],
    "step": f'{step_min}m'
  }

  if info['frame'] == 'GSE':
    data_jpl = get_horizons_coord(target_body, t).geocentricsolarecliptic.cartesian.xyz.to('km')
  if info['frame'] == 'GEI':
    data_jpl = get_horizons_coord(target_body, t).geocentricearthequatorial.cartesian.xyz.to('km')
  if info['frame'] == 'GSM':
    data_jpl = get_horizons_coord(target_body, t).geocentricsolarmagnetospheric.cartesian.xyz.to('km')

  resp = {
    'name': 'JPL',
    'time': time,
    'xyz': data_jpl.value.T/CFG['R_E'],
    'frame': info['frame']
  }

  with open(cache_file, 'wb') as f:
    print(f"Saving JPL data to cache file: {cache_file}")
    pickle.dump(resp, f)

  return resp


def _sscweb(info, logging=False):

  """
  1. SSCWeb provides TOD, J2K, GEO, GM, GSE, and SM
  TOD in SSCWeb HAPI server is same as GEI from SSCWeb.
  (The reason the HAPI server uses TOD for GEI is that "TOD" is the 
  POST query parameter names used to request data in GEI is "TOD"
  and consistency in the query parameters was sought.)

  2. From https://sscweb.gsfc.nasa.gov/users_guide/Appendix_C.html,

  > "Geocentric Equatorial Inertial system. This system has X-axis
  pointing from the Earth toward the first point of Aries (the position of
  the Sun at the vernal equinox). This direction is the intersection of the
  Earth's equatorial plane and the ecliptic plane and thus the X-axis lies
  in both planes. The Z-axis is parallel to the rotation axis of the Earth,
  and y completes the right-handed orthogonal set (Y = Z * X). Geocentric
  Inertial (GCI) and Earth-Centered Inertial (ECI) are the same as GEI."

  Based on the above quote, we treat GCI from CDAWeb as the same as GEI from SSCWeb.
  """
  server     = 'http://hapi-server.org/servers/SSCWeb/hapi'
  dataset    = info['dataset']
  frame      = info['frame']
  if frame == 'GEI':
    frame = 'TOD' # See SSCWeb note 1. above.
  if frame == 'GCI':
    frame = 'TOD' # See SSCWeb note 2. above.
  start      = info['start']
  stop       = info['stop']
  parameters = f'X_{frame},Y_{frame},Z_{frame}'
  opts       = {'logging': logging, 'usecache': True, 'cachedir': f"{CFG['cache_dir']}/positions/hapi"}

  data, meta = hapi(server, dataset, parameters, start, stop, **opts)

  info['xyz'] = np.column_stack((data[f'X_{frame}'], data[f'Y_{frame}'], data[f'Z_{frame}']))
  # Convert from YYYY-DOY to YYYY-MM-DD date format
  info['time'] = hapitime2datetime(data['Time'])

  # Return not needed b/c info is modified in place. Keep for clarity.
  return info


def _cdaweb(info, logging=False):

  if info is None:
    return None

  server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
  dataset    = info['dataset']
  parameter  = info['parameter']
  start      = info['start']
  stop       = info['stop']
  opts       = {'logging': logging, 'usecache': True, 'cachedir': f"{CFG['cache_dir']}/positions/hapi"}

  data, meta = hapi(server, dataset, parameter, start, stop, **opts)

  xyz = data[parameter]
  info['xyz'] = xyz/CFG['R_E'] # Convert from km to R_E
  # The HAPI SSCWeb server provides positions as three separate parameters. Here
  # we combine parameters into a list. SSCWeb reports in R_E while CDAWeb in km.
  # Convert CDAWeb data to R_E.
  info['time'] = hapitime2datetime(data['Time'])

  # Return not needed b/c info is modified in place. Keep for clarity.
  return info


def _amda(info, logging=False):

  if info is None:
    return None

  server     = 'https://amda.irap.omp.eu/service/hapi'
  dataset    = info['dataset']
  parameter  = info['parameter']
  start      = info['start']
  stop       = info['stop']
  opts       = {'logging': logging, 'usecache': True, 'cachedir': f"{CFG['cache_dir']}/positions/hapi"}

  data, meta = hapi(server, dataset, parameter, start, stop, **opts)

  xyz = data[parameter]
  info['xyz'] = xyz*CFG['R_E_AMDA']/CFG['R_E'] # Convert from AMDA Re to SSCWeb R_E
  info['time'] = hapitime2datetime(data['Time'])

  # Return not needed b/c info is modified in place. Keep for clarity.
  return info


def _print_first_last_(sscweb_, cdaweb_, amda_, jpl_):

  frame = sscweb_['frame']
  for i in [0, -1]:

    xyz_sscweb = sscweb_['xyz'][i]
    if cdaweb_ is not None:
      xyz_cdaweb = cdaweb_['xyz'][i]
    if jpl_ is not None:
      xyz_jpl = jpl_['xyz'][i]
    if amda_ is not None:
      xyz_amda = amda_['xyz'][i]

    time_sscweb = sscweb_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if cdaweb_ is not None:
      time_cdaweb = cdaweb_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if jpl_ is not None:
      time_jpl = jpl_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
    if amda_ is not None:
      time_amda = amda_['time'][i].strftime('%Y-%m-%dT%H:%M:%S.%fZ')

    u = 'R_E'
    if i == 0:
      print(f'  First values in {u}:')
    else:
      print(f'  Last values in {u}:')

    print('            {0:8s} {1:8s} {2:8s}'.format(f' X_{frame}', f' Y_{frame}', f' Z_{frame}'))
    print('    SSCWeb: {0:<8.4f} {1:<8.4f} {2:<8.4f} {3:s}'.format(*xyz_sscweb, time_sscweb))
    if cdaweb_ is not None:
      print('    CDAWeb: {0:<8.4f} {1:<8.4f} {2:<8.4f} {3:s}'.format(*xyz_cdaweb, time_cdaweb))
    if jpl_ is not None:
      print('       JPL: {0:<8.4f} {1:<8.4f} {2:<8.4f} {3:s}'.format(*xyz_jpl, time_jpl))
    if amda_ is not None:
      print('      AMDA: {0:<8.4f} {1:<8.4f} {2:<8.4f} {3:s}'.format(*xyz_amda, time_amda))


def _other_compare(satellite, infos, popts):
  from positions_plot import _plot
  # Compare SSCWeb J2K with SSCWeb GEI
  for info in infos['sscweb']:
    if info is not None and info.get('frame') == 'J2K':
      sscweb_j2k = _sscweb(info, logging=CFG['hapi_logging'])
    if info is not None and info.get('frame') == 'GEI':
      sscweb_gei = _sscweb(info, logging=CFG['hapi_logging'])

  _plot(satellite, sscweb_j2k, sscweb_gei, popts)


if __name__ == "__main__":
  from positions_plot import plot_combos

  args = cli()

  if len(args.satellites) == 0:
    print("No satellites specified. See python positions.py --help for more info.")
    exit(0)

  results = {}
  for satellite in args.satellites:

    infos = _infos(satellite)

    if False:
      _other_compare()

    results[satellite] = {}
    for i in range(len(infos['sscweb'])): # Loop over frames

      frame = infos['sscweb'][i]['frame']

      print(f"\n{satellite} {frame}")

      # Adds xyz and time to info dict
      sscweb = _sscweb(infos['sscweb'][i], logging=CFG['hapi_logging'])
      cdaweb = _cdaweb(infos['cdaweb'][i], logging=CFG['hapi_logging'])
      amda = _amda(infos['amda'][i], logging=CFG['hapi_logging'])

      times = infos['sscweb'][i]['time']
      jpl = _jpl(infos['jpl'][i], times, logging=CFG['hapi_logging'])
      _print_first_last_(sscweb, cdaweb, amda, jpl)

      sources = [sscweb, cdaweb, amda, jpl]
      results[satellite][frame] = sources

      plot_combos(satellite, sources, CFG['plot_opts'])

  if args.satellites == CFG['known_satellites']:
    fname = f"{CFG['cache_dir']}/positions/positions"
  else:
    fname = f"{CFG['cache_dir']}/positions/positions-{'_'.join(args.satellites)}"

  print(f"Writing positions data to {fname}.{{json,pkl}}")
  utilrsw.write(f'{fname}.json', results)
  utilrsw.write(f'{fname}.pkl', results)
