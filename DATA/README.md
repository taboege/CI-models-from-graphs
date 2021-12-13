# Quick guide to the files

## Simplicial complexes in `Macaulay2` notation, one per line

Arbitrary simplicial complexes:

- complexes-3.txt
- complexes-4.txt
- complexes-5.txt
- complexes-6.txt

Connected matroids (the ones on 8 vertices are up to rank 4):

- matroids-4c.txt
- matroids-5c.txt
- matroids-6c.txt
- matroids-7c.txt
- matroids-8c-r4.txt

**Caution:** currently the ground set for matroids is indexed starting
from **0** while general complexes start from **1**. To be fixed soon.

## Dimension, degree and Betti table

Corresponding to the inputs above, there are listings of triples
`(dim, codim, deg, betti)` of the dimension, codimension, degree
and the entire Betti table as printed by `Macaulay2` for the
toric variety:

- complexes-3-betti-mingens.txt
- complexes-4-betti-mingens.txt
- complexes-5-betti-mingens.txt
- complexes-6-betti-mingens.txt
- matroids-4c-betti-mingens.txt
- matroids-5c-betti-mingens.txt
- matroids-6c-betti-mingens.txt
- matroids-7c-betti-mingens.txt
- matroids-8c-r4-betti-mingens.txt

## aED, pED and ML degrees

**These numbers still need to be confirmed:**

Euclidean distance degrees, affine and projective.
Each pair `affine(c,r)` gives the degree `c` and the number of
real solutions `r` for **some** random sample.

- complexes-3-EDdeg.txt
- complexes-4-EDdeg.txt
- complexes-5-EDdeg.txt

Maximum likelihood degree. They are written in the same format
as above but there is only one sort of degree (which is always
called "affine"):

- complexes-3-MLdeg.txt
- complexes-4-MLdeg.txt
- complexes-5-MLdeg.txt
