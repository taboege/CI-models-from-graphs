============================
Marginal Independence Models
============================

| This page contains selected examples and the database of marginal independence models from the paper:
| Tobias Boege, Sonja Petrović, Bernd Sturmfels: Marginal Independence Models.
.. | ARXIV: https://arxiv.org/abs/...
| CODE/DATA: https://mathrepo.mis.mpg.de/MarginalIndependence/

ABSTRACT: We impose rank one constraints on marginalizations of a tensor, given by a simplicial complex. Following work of Kirkup and Sullivant, such marginal independence models can be made toric by a linear change of coordinates. We study their toric ideals, with emphasis on random graph models and independent set polytopes of matroids. We develop the numerical algebra of parameter estimation, using both Euclidean distance and maximum likelihood, and we present a comprehensive database of small models.

We consider :math:`n` discrete random variables :math:`X_1, X_2, \ldots, X_n` where :math:`X_i` takes
values in the finite set :math:`[d_i] = \{ 1, 2, \ldots, d_i \}`. A joint distribution for these random
variables is a tensor :math:`P = (p_{i_1 i_2 \cdots i_n})` of format :math:`d_1 \times d_2 \times
\cdots \times d_n` whose entries are nonnegative real numbers that sum to :math:`1`.

Every subset :math:`\sigma \subseteq [n]` defines a *marginal* distribution :math:`P_\sigma` by summing
out all dimensions of the tensor outside of :math:`\sigma`. This is a tensor of format :math:`\bigtimes_{i \in \sigma} d_i`.
The random variables indexed by :math:`\sigma` are *completely independent* if the marginal :math:`P_\sigma`
has rank :math:`1`. This notion of complete independence is hereditary: if :math:`P_\sigma` has rank :math:`1`,
then all subsets of :math:`\sigma` have the same property. Hence the set of completely independent
subcollections of any tensor :math:`P` is a simplicial complex.

Let, conversely, :math:`\Sigma` be any simplicial complex with vertex set :math:`[n]`. We assume
:math:`\{i\} \in \Sigma` for all :math:`i \in [n]`. The *marginal independence model* :math:`\mathcal{M}_\Sigma`
is the set of distributions :math:`P \in \Delta_{D-1}` such that the marginal tensor :math:`P_\sigma`
has rank :math:`1` for all :math:`\sigma \in \Sigma`.

The rank :math:`1` constraints translate to vanishing of :math:`2\times 2` minors of flattenings of the
tensor :math:`P`. This gives the ideal of the marginal independence model in *natural joint probability
coordinates*, i.e., the entries of :math:`P`. Its generators are quadratic polynomials.
Under a linear change of coordinates, the ideal becomes toric, as proved in the paper. The new
coordinates are called *Möbius coordinates*. This is implemented in our ``Macaulay2`` package
``MarginalIndependenceModels``.

The following example demonstrates the usage of this package at the binary :math:`3`-cycle,
i.e., the simplicial complex with facets :math:`\{ 12, 13, 23 \}` imposing independence
conditions on three binary random variables. The example also appears in S. Sullivant:
*Algebraic Statistics* (2018) and G. Kirkup: *Random variables with completely independent
subcollections* (2007).

.. code-block:: macaulay2

                needsPackage "MarginalIndependenceModels"
                sigma = {{1,2},{1,3},{2,3}}
                I = marginalIndependenceIdeal sigma

The matrix defining this toric ideal is given by the following polytope:

.. code-block:: macaulay2

                P = binaryPolytope sigma
                vertices P

To appreciate the linear change  of coordinates, take a look at the map:

.. code-block:: macaulay2

                phi = moebius sigma

The target and source rings are in natural probabilities and relevant Moebius coordinates, respectively.
***SAY A FEW MORE WORDS IN DETAIL?...***

.. code-block:: macaulay2

                target phi
                source phi

Caveat
It is worth noting that in this package, the ambient polynomial ring uses only the relevant
Möbius coordinates, while the full model lives in a larger space. Hence there will be a
discrepancy comparing the theoretical model results and properties of the ideal :math:`I`,
for example, dimension is offset by the number of irrelevant Möbius coordinates.



TODO: (1) aED, pED, ML degree computations using ``HomotopyContinuation.jl``.
(2) Include Bernd's code for Example 1 (the ternary 3-cycle Kirkup couldn't compute).
(39 Probably show how to do one challenging examples with the ``betti mingens``,
like Example 23 (real projective plane).



Project page created: 14/12/2021

Project contributors: Tobias Boege, Sonja Petrović, Bernd Sturmfels

Corresponding author of this page: Tobias Boege, tobias.boege@mis.mpg.de
