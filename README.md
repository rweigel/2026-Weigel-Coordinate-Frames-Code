[![DOI](https://zenodo.org/badge/1241504388.svg)](https://doi.org/10.5281/zenodo.20347420)

Software used for "Coordinate Systems and Transforms in Space Physics: Terms, Definitions, Implementations, and Recommendations for Reproducibility by
R.S. Weigel, A.Y. Shih, R. Ringuette, I. Christopher, S.M. Petrinec, S. Turner, R.M. Candey, G.K. Stephens, B. Cecconi (https://arxiv.org/abs/2601.07605).

To install, use

```
git clone https://github.com/rweigel/2026-Weigel-Coordinate-Frames-Code
git clone https://github.com/rweigel/2026-Weigel-Coordinate-Frames-Data

cd 2026-Weigel-Coordinate-Frames-Code
pip install -e .
```

To generate the figures in the paper using cached data in `2026-Weigel-Coordinate-Frames-Data`, use

```
python angles_plot.py
python positions_plot.py --satellites Geotail,MMS-2
```

To re-download data and perform computations for section 4.1, remove or rename the directory `2026-Weigel-Coordinate-Frames-Data/positions` and execute

```
python positions.py --satellites Geotail,MMS-2
```

To compute angles for section 4.2, remove or rename the directory `2026-Weigel-Coordinate-Frames-Data/angles` and execute

```
python angles.py
```
