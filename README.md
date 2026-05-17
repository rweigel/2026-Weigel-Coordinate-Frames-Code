```
git clone https://github.com/rweigel/2026-Weigel-Coordinate-Frames-Code
git clone https://github.com/rweigel/2026-Weigel-Coordinate-Frames-Data

cd 2026-Weigel-Coordinate-Frames-Code
pip install -e .
```

To generate figures in paper using cached data in `2026-Weigel-Coordinate-Frames-Data`, use

```
python angles_plot.py
python positions_plot.py --satellites Geotail,MMS-2
``

To re-download data and perform computations for section 4.1, remove or rename `2026-Weigel-Coordinate-Frames-Data/positions` and execute
```
python positions.py --satellites Geotail,MMS-2
```

To compute angles for section 4.2, remove or rename `2026-Weigel-Coordinate-Frames-Data/angles` and execute
```
python angles.py
```
