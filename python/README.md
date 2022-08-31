About
=====

This is a simplified Numpy implementation that implements a tree hierarchy, but
uses a per-leaf list of remote boxes to determine how to accumulate remote
vectors, for simplicity.  The cos and sin vectors are generated recursively
using the addition formulas.

Prerequisites
=============

Numpy and Scipy are required in order for this to work; both are readily
availble through [Conda](https://docs.conda.io/en/latest/) or
[PIP](https://pypi.org/project/pip/).  The folder of spherical t-designs should
first be installed, as detailed on this project's homepage.
