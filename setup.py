import pathlib
import setuptools

install_requires = [
  'numpy',
  'pandas',
  'matplotlib',
  'hapiclient',
  'sunpy>=6.1.0', # https://github.com/sunpy/sunpy/pull/8193
  'datetick @ git+https://github.com/rweigel/datetick@main',
  'hxform @ git+https://github.com/rweigel/hxform@main',
  'utilrsw[time] @ git+https://github.com/rweigel/utilrsw@main',
  'utilrsw[mpl] @ git+https://github.com/rweigel/utilrsw@main'
]


readme = pathlib.Path(__file__).with_name('README.md')
long_description = readme.read_text(encoding='utf-8') if readme.exists() else ''

kwargs = {
  'name': '2026-Weigel-Coordinate-Frames-Code',
  'author': 'Bob Weigel',
  'author_email': 'rweigel@gmu.edu',
  'packages': setuptools.find_packages(),
  'description': 'Code for Weigel et al. 2026, https://arxiv.org/abs/2401.07605',
  'long_description': long_description,
  'long_description_content_type': 'text/markdown',
  'include_package_data': True,
  'package_data': {'': ['README.md']},
  'license': 'LICENSE.txt',
  'install_requires': install_requires
}

setuptools.setup(**kwargs)
