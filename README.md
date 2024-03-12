# Small World Propensity
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10299681.svg)](https://doi.org/10.5281/zenodo.10299681)
[![GitHub release](https://img.shields.io/github/v/release/rkdan/small_world_propensity?include_prereleases)](https://GitHub.com/rkdan/small_world_propensity/releases)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/small-world-propensity.svg)](https://pypi.python.org/pypi/small-world-propensity/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

This python package was adapted from the MATLAB package as first presented in [Small-World Propensity and Weighted Brain Networks](https://www.nature.com/articles/srep22057) (2016) by Sarah Feldt Muldoon, Eric W. Bridgeford & Danielle S. Bassett. Their original MATLAB implementation can be found [here](https://kk1995.github.io/BauerLab/BauerLab/MATLAB/lib/+mouse/+graph/smallWorldPropensity.html).

## Use
The small-world propensity package can be installed using pip
```
python -m pip install small-world-propensity
```
`small_world_propensity` can be called in two ways: either with a single adjacency matrix, or with a list of adjacency matrices and a boolean list denoting whether each matrix is binary or not. In either case, `small_world_propensity` will return a `pandas` dataframe similar to the following:
![Dataframe](https://github.com/rkdan/small_world_propensity/blob/main/img/dataframe.png?raw=True)

## Generation of regular and random matrices
Using the structural network of the cat cortex obtained from tract-tracing studies between 52 brain regions, we can visualize the process behind the calculation of the small-world propensity, $\phi$. The matrix is loaded using

```
cat = sio.loadmat('data/cat.mat')['CIJctx']
```
We can then ensure symmetry by calling
```
symm_cat = swp.make_symmetric(cat)
```
In order to get the regular version of the cat matrix, we first find the effective average radius:
```
r = swp.get_avg_rad_eff(symm_cat)
cat_reg = swp.regular_matrix_generator(symm_cat, r)
```
Finally we produce the randomized cat matrix:
```
cat_rand = swp.randomize_matrix(cat_symm)
```
![Cat matrices](https://github.com/rkdan/small_world_propensity/blob/main/img/cat.png?raw=True)

The graphs visualized in a circular layout look as follows:

![Cat graphs](https://github.com/rkdan/small_world_propensity/blob/main/img/cat_graphs.png?raw=True)

## Comparison of $\phi$ in real networks
We can take the networks used in _Muldoon et al_ and plot $\phi$, $\Delta_L$, $\Delta_C$, and $\delta$. Note that these networks are not the exact same as the ones used in _Muldoon et al_, and due to differences in how Numpy performs permutations, and the use of NetworkX and iGraph libraries, the results are not identical, but still match closely.

The adjacency matrices:
![Adjacency matrices](https://github.com/rkdan/small_world_propensity/blob/main/img/matrices.png?raw=True)

And the results:
![Summary](https://github.com/rkdan/small_world_propensity/blob/main/img/summary.png?raw=True)

To cite this work, please use:
```bibtex
@software{small-world-propensity,
  author       = {{Daniels, R. K.}},
  title        = {small-world-propensity},
  year         = 2023,
  publisher    = {Zenodo},
  version      = {v0.0.8},
  doi          = {10.5281/zenodo.10299681},
  url          = {https://github.com/rkdan/small-world-propensity}
}
```
Please also cite the authors of the original MATLAB implementation:
```bibtex
@article{Muldoon2016,
    author = "Muldoon, Sarah Feldt and Bridgeford, Eric W. and Bassett, Danielle S.",
    title = "{Small-World Propensity and Weighted Brain Networks}",
    doi = "10.1038/srep22057",
    journal = "Scientific Reports",
    volume = "6",
    number = "1",
    pages = "P07027",
    year = "2016"
}
```

> [!NOTE]  
> This software has a GNU AGPL license. If this license is inadequate for your use, please get in touch.