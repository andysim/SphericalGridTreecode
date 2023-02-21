C++ SGT implementation example
==============================

This is a simple C++ example that uses and upward and a downward pass to
accumulate / distribute the structure factors.  The calculation is parallelized
by distributing over the quadrature points for simplicity, with the resulting
energies and forces from each core accumulated at the end.

Prerequisites
-------------

As noted in the parent directory's notes, we assume that the `SS31-Mar-2016`
folder has been downloaded and placed one level up from this directory.

Building
--------

This example code uses CMake to build, and will automatically detect OpenMP for
parallelism.  This code uses special functions from a header-only subset of
Boost C++ library, which is bundled as part of this codebase; CMake will
automatically set the build up to find these files.  To build, run:

``` bash
mkdir build
cd build
cmake ..
make -j8
```

If the build ran successfully then, _e.g._, to run a waterbox using a real
space cutoff of 12.0, aiming for 1E-9 precision on 4 cores simply run:

``` bash
./run_waterbox 9 12.0 4 ../../SS31-Mar-2016 ../Data/waterbox.xyz
```

The waterbox xyz file follows the Tinker format: another simple utility is
provided to build larger waterboxes by translating along the x, y, and/or, z
axes.
