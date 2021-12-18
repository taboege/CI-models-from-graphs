============================
Marginal Independence Models
============================

| This page contains selected examples and the database of marginal independence models from the paper:
| Tobias Boege, Sonja Petrović, Bernd Sturmfels: Marginal Independence Models.
.. | ARXIV: https://arxiv.org/abs/...
| CODE: https://github.com/taboege/CI-models-from-graphs/blob/master/MarginalIndependenceModels.m2
| DATA: https://github.com/taboege/CI-models-from-graphs/tree/master/DATA
| SEE ALSO: https://mathrepo.mis.mpg.de/MarginalIndependence/

ABSTRACT: We impose rank one constraints on marginalizations of a tensor, given by a simplicial complex.
Following work of Kirkup and Sullivant, such marginal independence models can be made toric by a linear
change of coordinates. We study their toric ideals, with emphasis on random graph models and independent
set polytopes of matroids. We develop the numerical algebra of parameter estimation, using both Euclidean
distance and maximum likelihood, and we present a comprehensive database of small models.

We consider :math:`n` discrete random variables :math:`X_1, X_2, \ldots, X_n` where :math:`X_i` takes
values in the finite set :math:`[d_i] = \{ 1, 2, \ldots, d_i \}`. A joint distribution for these random
variables is a tensor :math:`P = (p_{i_1 i_2 \cdots i_n})` of format :math:`d_1 \times d_2 \times
\cdots \times d_n` whose entries are nonnegative real numbers that sum to :math:`1`.

Every subset :math:`\sigma \subseteq [n]` defines a *marginal* distribution :math:`P_\sigma` by summing
out all dimensions of the tensor outside of :math:`\sigma`. This is a tensor of format :math:`\bigtimes_{i \in \sigma} d_i`
which represents the distribution of the subset of variables :math:`X_i`, :math:`i \in \sigma`.
The random variables indexed by :math:`\sigma` are *completely independent* if the marginal tensor
:math:`P_\sigma` has rank :math:`1`. Rank :math:`0` is impossible due to the sum to one condition,
so marginal independence is equivalent to the tensor having the minimal rank possible. This notion
of complete independence is hereditary: if :math:`P_\sigma` has rank :math:`1`, then all subsets of
:math:`\sigma` have the same property. Hence the set of completely independent subcollections of
any tensor :math:`P` is a simplicial complex. These complexes always include the singletons
:math:`\{i\}` whose marginal tensors are just non-zero vectors.

Let, conversely, :math:`\Sigma` be any simplicial complex with vertex set :math:`[n]`. We assume
:math:`\{i\} \in \Sigma` for all :math:`i \in [n]`. The *marginal independence model* :math:`\mathcal{M}_\Sigma`
is the set of distributions :math:`P \in \Delta_{D-1}` such that the marginal tensor :math:`P_\sigma`
has rank :math:`1` for all :math:`\sigma \in \Sigma`.

The rank :math:`1` constraints translate to vanishing of all :math:`2 \times 2` minors of flattenings
of each marginal  :math:`P_\sigma`. This gives the ideal of the marginal independence model in *natural
joint probability coordinates*, i.e., the entries of :math:`P`. Its generators are quadratic polynomials.
Under a linear change of coordinates, the ideal becomes toric, as proved in the paper. The new
coordinates are called *Möbius coordinates*. This is implemented in our Macaulay2_ package
MarginalIndependenceModels_.

Polytope, ideal, dimension and degree
-------------------------------------

The following example demonstrates the usage of this package at the binary :math:`3`-cycle,
i.e., the simplicial complex with facets :math:`\{ 12, 13, 23 \}` imposing independence
conditions on three binary random variables. The example also appears in S. Sullivant:
*Algebraic Statistics* (2018) and G. Kirkup: *Random variables with completely independent
subcollections* (2007).

.. code-block:: macaulay2

                needsPackage "MarginalIndependenceModels"
                Sigma = {{1,2},{1,3},{2,3}}
                I = marginalIndependenceIdeal sigma

It is enough to specify the simplicial complex ``Sigma`` in terms of its facets;
the set of all faces is computed on the fly if needed. The matrix defining this
toric ideal is given by the following polytope:

.. code-block:: macaulay2

                needsPackage "Polyhedra"
                P = binaryPolytope Sigma
                vertices P

To appreciate the linear change of coordinates, take a look at the map:

.. code-block:: macaulay2

                phi = moebius Sigma

The target and source rings are in natural probabilities and relevant Möbius
coordinates, respectively.

.. code-block:: macaulay2

                target phi
                source phi

The marginal tensor :math:`P_{12}` is a :math:`2 \times 2` matrix which must have rank
:math:`1` in the above 3-cycle model. This is imposed in Möbius coordinates by the
quadric :math:`q q_{12} - q_1 q_2`. In tensor coordinates, this is the quadratic
non-binomial :math:`a`, as confirmed by

.. code-block:: macaulay2

                preimage_phi( q_{}*q_{1,2} - q_{1}*q_{2} )

Essential geometric information about this model can be computed and checked against
the Kirkup's paper:

.. code-block:: macaulay2

                codim I, degree I, betti mingens I

**Caveat:** The ambient polynomial ring of ``I`` uses only the relevant Möbius coordinates,
while the full model lives in a larger space. Hence there will be a discrepancy when computing
the model dimension in Macaulay2. The dimension of the model as a projective variety can
be obtained as ``2^3 - 1 - codim I``.

We refer to the documentation of the package for more information.

Minimal generators for the real projective plane model
------------------------------------------------------

Example 23 is now quick and easy to confirm:

.. code-block:: macaulay2

                RP = {{1,2,3},{1,2,4},{1,3,5},{1,4,6},{1,5,6},{2,3,6},{2,4,5},{2,5,6},{3,4,5},{3,4,6}}
                betti mingens marginalIndependenceIdeal RP

The ternary 3-cycle
-------------------

Our code only covers models on binary random variables. Checking the free resolution of
the ternary model of the 3-cycle in Example 1 requires manual coding:

.. code-block:: macaulay2

                R = QQ[
                  -- Parameters for the Segre variety
                  z, a1, b1, c1, a2, b2, c2,
                  -- Relevant Moebius coordinates, S stands for "sum"
                  q11S, q12S, q21S, q22S, q1S1, q1S2, q2S1, q2S2,
                  qS11, qS12, qS21, qS22,
                  q1SS, q2SS, qS1S, qS2S, qSS1, qSS2, qSSS,
                  MonomialOrder => Eliminate 7
                ];

                -- Parametrization of the relevant Moebius coordinates as in Section 3.
                I = ideal(
                  q11S - z*a1*b1*(1-c1-c2),
                  q12S - z*a1*b2*(1-c1-c2),
                  q21S - z*a2*b1*(1-c1-c2),
                  q22S - z*a2*b2*(1-c1-c2),
                  -- etc.
                );

                -- Implicitize
                M = eliminate({z, a1, b1, c1, a2, b2, c2}, I);

                betti mingens M
                betti res M

The entire file is attached to this page: Cycle3-3x3x3.m2_.

Euclidean distance and maximum likelihood degrees
-------------------------------------------------

To compute aED, pED and ML degrees of our models, we used the Julia_ package HomotopyContinuation.jl_.
We employ the parametrization of the model in Möbius coordinates to write the critical equations
for the unconstrained optimization problem. In these coordinates, marginal independences correspond
to factorizations of the relevant Möbius coordinates into the singletons:
:math:`q_\sigma = \prod_{i \in \sigma} q_i`.

.. code-block:: julia

                using HomotopyContinuation

                @var q1 q2 q3 q12 q13 q23 q123 z

                sample = randn(ComplexF64, 8)
                p000,p001,p010,p011,
                p100,p101,p110,p111 = sample

                diffs = [
                  -p000 + z*(q123),
                  -p001 + z*(q12-q123),
                  -p010 + z*(q13-q123),
                  -p100 + z*(q23-q123),
                  -p101 + z*(q2-q12-q23+q123),
                  -p110 + z*(q3-q13-q23+q123),
                  -p011 + z*(q1-q12-q13+q123),
                  -p111 + z*(1-q1-q2+q12-q3+q13+q23-q123)
                ]

                dist = sum([d^2 for d in diffs])

This sets up the parametrization and the objective function.
The 3-cycle corresponds to the below factorization:

.. code-block:: julia

                subslist = [ z=>1, q12=>q1*q2, q13=>q1*q3, q23=>q2*q3 ]
                dist = subs(dist, subslist...)

HomotopyContinuation.jl can write out the critical equations and compute the number of complex
solutions for the random complex data chosen above.

.. code-block:: julia

                vars = variables(dist)
                eqns = differentiate(dist, vars)

                R = solve(eqns; show_progress = false,
                  tracker_options = TrackerOptions(
                    automatic_differentiation = 3,
                    parameters = :conservative
                  )
                )

                C = certify(eqns, R; show_progress = false)
                println(ndistinct_certified(C))

To change the model, only the variable ``subslist`` has to be altered. To automate the
computation of degrees, the above program was embedded into a Perl_ program which
modifies the ``subslist`` for a given simplicial complex, runs the Julia program and
parses its output. The parametrizations are obtained via our Macaulay2 code and hardcoded
in each file:

- <EDdeg3.jl.pl>
- <EDdeg4.jl.pl>
- <EDdeg5.jl.pl>
- <MLdeg3.jl.pl>
- <MLdeg4.jl.pl>

Plain text database of complexes, models and degrees
----------------------------------------------------

Our findings are summarized in plain text database files. They list each simplicial complex
and the data for its model in one line. This includes the facet description of the complex,
its dimension, number of faces, facets and f-vector to make it easy to find a given complex
by hand. We indicate whether the complex is a matroid or even the complex of forests of an
undirected simple graph. For its toric ideal we list the dimension, codimension and degree,
moreover its m-vector, i.e., for each degree the number of minimal generators. Finally,
whenever available, we give the aED, pED and ML degree. The complexes are listed up to
isomorphy and they always contain :math:`\{i\}` for all :math:`i \in [n]`.

- <MarginalIndependence-3.txt>
- <MarginalIndependence-4.txt>

.. - <MarginalIndependence-5.txt>
.. 
.. For higher :math:`n`, we are only able to list combinatorial and ideal data for models
.. induced by matroids:
.. 
.. - <MatroidIndependence-6.txt>
.. - <MatroidIndependence-7.txt>

Colophon
--------

Project page created: 18/12/2021

Project contributors: Tobias Boege, Sonja Petrović, Bernd Sturmfels

Corresponding author of this page: Tobias Boege, tobias.boege@mis.mpg.de

.. _Macaulay: https://faculty.math.illinois.edu/Macaulay2/
.. _MarginalIndependenceModels: https://github.com/taboege/CI-models-from-graphs/blob/master/MarginalIndependenceModels.m2
.. _Julia: https://julialang.org
.. _HomotopyContinuation.jl: https://www.juliahomotopycontinuation.org
.. _Perl: https://www.perl.org

.. _Cycle3-3x3x3.m2: Cycle3-3x3x3.m2
