import os

CFG = {
  # Usedy by angles*py and positions*.py
  'cache_dir': os.path.join(os.path.dirname(__file__), '..', '2026-Weigel-Coordinate-Frames-Data'),

  # The following are used only by positions*.py
  'hapi_logging': False,

  # Earth radius in km used by SSCWeb. See
  # https://sscweb.gsfc.nasa.gov/users_guide/Users_Guide_pt1.html#1.3
  # Source code for SSCWeb is in private STCT Google Drive folder in file
  # icss_rel82.tar.gz (ISTP coordinate transformation library). It is not, and
  # cannot be, publicly available).
  'R_E': 6378.16, # km

  # In https://amda.irap.omp.eu/service/hapi/info?id=clust1-orb-all
  # the referenced SPASE record says "Units": "Re", "UnitsConversion": "6.4e6>m"
  # From Vincent Génot to Weigel on 2025-05-15: SPASE is wrong; correct value is:
  'R_E_AMDA': 6378.137, # km

  'known_satellites': [
    'Cluster-1',
    'Cluster-2',
    'Cluster-3',
    'Cluster-4',
    'DSCOVR',
    'Geotail',
    'THEMIS-A',
    'THEMIS-B',
    'THEMIS-C',
    'THEMIS-D',
    'MMS-1',
    'MMS-2',
    'MMS-3',
    'MMS-4'
  ],
}

# Options for positions_plot.py
CFG['plot_opts'] = {
  'R_E': CFG['R_E'],
  'cos_warnings': False,

  # Interpolate to union of timestamps from both sources before comparing.
  # If False, then only compare values at common timestamps.
  'interp_times': True,
  # If interp_times is True, then only compare values at timestamps that are within this tolerance of each other.
  'interp_times_tolerance': 10, # seconds

  'print_interp_vals': False, # Print interpolated values
  'print_common_times': False, # Print vals for common times

  'fontsize': 20,
  'figsize_inches': (170/25.4, 228/25.4), # (width, height) in inches
  'legend_kwargs': {
    'columnspacing': 0.65,
    'handletextpad': 0.5,
    'handlelength': 1.5,
    'borderaxespad': 0.0,
    'loc': 'upper center'
  }

}


def cli():
  import argparse

  def parse_satellites(value):
    return [part.strip() for part in value.split(',') if part.strip()]

  # Check that provided satellites are valid
  def validate_satellites(satellites):
    for satellite in satellites:
      if satellite not in CFG['known_satellites']:
        raise argparse.ArgumentTypeError(f"Invalid satellite: {satellite}. Must be one of {CFG['known_satellites']}")
    return satellites

  description="Compute and plot magnetic field models for satellite data."
  parser = argparse.ArgumentParser(description=description)

  parser.add_argument(
    "--satellites",
    dest="satellites",
    type=parse_satellites,
    help=f"Comma-separated list of external satellites. Possible values {CFG['known_satellites']}. Default is all satellites.",
    default=CFG['known_satellites']
  )

  args = parser.parse_args()
  validate_satellites(args.satellites)

  return args

