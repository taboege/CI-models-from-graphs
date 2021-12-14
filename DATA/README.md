# Quick guide to the files

## Simplicial complexes in `Macaulay2` notation, one per line

Arbitrary simplicial complexes:

- complexes-3.txt
- complexes-4.txt
- complexes-5.txt
- complexes-6.txt

Matroids (loopless, up to isomorphy) go up to seven ground elements:

- matroids-4.txt
- matroids-5.txt
- matroids-6.txt
- matroids-7.txt

**Caution:** it is possible that the isomorphy representatives chosen
for matroids are not the same as the ones chosen for complexes!

## Dimension, degree and Betti table

Corresponding to the inputs above, there are listings of triples
`(dim, codim, deg, betti)` of the dimension, codimension, degree
and the entire Betti table as printed by `Macaulay2` for the
toric variety:

- complexes-3-betti-mingens.txt
- complexes-4-betti-mingens.txt
- complexes-5-betti-mingens.txt
- complexes-6-betti-mingens.txt
- matroids-4-betti-mingens.txt
- matroids-5-betti-mingens.txt
- matroids-6-betti-mingens.txt
- matroids-7-betti-mingens.txt

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
