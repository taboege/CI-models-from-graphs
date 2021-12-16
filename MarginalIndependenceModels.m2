newPackage("MarginalIndependenceModels",
    Version => "1.0",
    Date => "December 13, 2021",
    Authors => {
	{
	    Name => "Tobias Boege",
	    Email => "post@taboege.de",
	    HomePage => "taboege.de"
	},
        {
	    Name => "Contributing Author: Sonja Petrovic",
	    Email => "sonja.petrovic@iit.edu",
	    HomePage => "sonjapetrovicstats.com"
	}
    },
    Headline => "Toric ideals of marginal independence models for n-way tensors"
  )


export {
    "binaryPolytope",
    "moebius",
    -- "quadraticModel", 
    "marginalIndependenceIdeal"
    }

needsPackage "FourTiTwo"; --we compute a Markov basis
needsPackage "Polyhedra"; --we construct a polytope

-- the IndexedVariable names for two coordinate rings: 
p = local p;
q = local q;



--***************************************--
--  Exported methods 	     	     	 --
--***************************************--

-------------------------------------
--      Model polytope, binary case 
-------------------------------------
-*
Input:
  FF = a simplicial complex, given as a list of facets
Output:
  The 0/1 polytope whose vertices are the faces of FF, given as a polyhedron.
  (Note: this is a type defined in "Polyhedra" package.)
*-
binaryPolytope = method(
    TypicalValue => Polyhedron
    );
binaryPolytope (List) := Polyhedron => FF -> (
  N := sort unique flatten FF;
  FF = complexify FF; -- really important in this case
  convexHull transpose matrix(
    apply(
      select(subsets(N), A -> member(A, FF)),
      -- This adds a row of homogenizing 1s at the end
      A -> flatten { apply(N, i -> if member(i, A) then 1 else 0), 1 }
    )
  )
);




----------------------------------------
--      Moebius setup
----------------------------------------
-*
Input:
  FF = simplicial complex, given as a list of facets
Output:
  The linear ring map which transforms the natural probability variables p into relevant Moebius coordinates q. 
*-
moebius = method(
    TypicalValue => RingMap
    );
moebius (List) := RingMap => FF -> (
  q = symbol q;
  p = symbol p;
  N := sort unique flatten FF;
  Q := QQ[ apply(moebiusLabels(N), A -> q_A) ];
  P := QQ[ apply(binaryLabels(N),  a -> p_a) ];
  map(P, Q, apply(moebiusLabels(N), A -> q_A => (
    -- D is the set where 0/1 should vary in the summation
    D := select(N, d -> not member(d, A));
    sum apply(apply(subsets(D), B -> (
      a := new MutableList from #N:0;
      apply(D, d -> (
        a#(position(N, i -> i == d)) = if member(d, B) then 1 else 0;
      ));
      toSequence a
    )), a -> p_a)
  )))
);


--------------------------------------------------------------------------------
--      The ideal of the model of marginal independence
--------------------------------------------------------------------------------
-*
Input:
  FF = simplicial complex, given as a list of facets
Output:
  a toric ideal of the model of marginal independence in Moebius coordinates 
*-
marginalIndependenceIdeal = method(
    TypicalValue => Ideal --ideal or list? or matrix of gens? 
    );    
marginalIndependenceIdeal (List) := Ideal => FF -> (
  q = symbol q; -- this is guided by how indexed variables were handled in GraphicalModels.m2 :) 
  N := sort unique flatten FF;
  B := binaryPolytope FF;
  FF = complexify FF;
  Q := QQ[ select(subsets(N), A -> member(A, FF)) / sort / (A -> q_A) ];
  toricMarkov(Polyhedra$vertices(B), Q)
);



-----------------------------------------------------
--      Some distinguished quadrics in the model
-----------------------------------------------------
-*
Input:
  FF = simplicial complex, given as a list of facets
Output:
  A subideal of the marginal independence ideal generated by quadrics of the form ****WRITE ME CONCISELY****. 
  This contains only some special quadrics that vanish on the model, but may be useful to some.
*-
quadraticModel = method(
    TypicalValue => Ideal
    );    
quadraticModel (List) := Ideal => FF -> (
  q = symbol q;    
  N := sort unique flatten FF;
  f := moebius(FF);
  Q := source f; use Q;
  FF = complexify(FF) / set;
  good := (s,t) -> member(set(s) + set(t), FF);
  ideal flatten apply(subsets(N), s -> (
    select(unique apply(select(subsets(N), t -> good(s,t)), t -> (
      i := sort toList(set(s) * set(t));
      u := sort toList(set(s) + set(t));
      q_s * q_t - q_i * q_u
    )), f -> f != 0)
  ))
);


--***************************************************************************--
--                          Internal functions 	                     	     	 --
--***************************************************************************--

--------------------------------------------------------------------------------
--      Simplicial complex and independence list housekeeping functions
--------------------------------------------------------------------------------

-- Generate the simplicial complex from a sequence of sets
complexify = method()
complexify (List) := List => FF -> unique flatten apply(FF, F -> subsets(F) / sort)

-- Conditional independence statements corresponding to the
-- independences in the complex FF.
CIstmts = N -> (
  flatten apply(subsets(N, 2) / sort, ij -> (
    M := toList(difference(set N, ij));
    apply(subsets(M) / sort, K -> {ij#0, ij#1, K})
  ))
);
indepCI = FF -> unique flatten apply(FF, CIstmts); 



--------------------------------------------------------------------------------
--      Labels for two rings to create the Moebius coordinate change map
--------------------------------------------------------------------------------
-*
Input:
  N = simplicial complex given as a list on n elements
Output:
  (binary labels): list indices of natural joint probabilities on the n random variables. 
    These are used as indices of the IndexedVariable "p".
  or 
  (moebius labels): list of indices of relevant Moebius coordinates for the marginal independence model given by N. 
    These are used as indices of IndexedVariable "q". 
*-
binaryLabels = method(
    TypicalValue => List
    );
binaryLabels (List) := List => N -> (
  A := set{0,1};
  toList(fold((X,Y) -> (X ** Y) / splice, #N:A)) / toSequence
);

moebiusLabels = method(
    TypicalValue => List
    );
moebiusLabels (List) := List => N -> (
  subsets(N)
);




--***************************************************************************--
--                          DOCUMENTATION     	                     	     	 --
--***************************************************************************--
beginDocumentation()
doc ///
  Key
    MarginalIndependenceModels
  Headline
    a package for computing marginal independence models for n-way tensors
  Description
    Text
      Let $X_1,X_2,\ldots,X_n$ be random variables, where $X_i$ has the finite state space $[d_i] = \{1,2,\ldots,d_i \}$.
      A joint distribution for these random variables is a tensor $P = (p_{i_1 i_2 \cdots i_n})$ 
      of format $d_1 \times d_2 \times \cdots  \times d_n$ whose entries are nonnegative real numbers that sum to $1$. 

      Let $\Sigma$ be any simplicial complex with vertex set $[n] $. We assume $\{i\} \in \Sigma$ for all $i \in [n]$.
      The {\em marginal independence model} $\mathcal{M}_\Sigma$  is the set of distributions $P \in \Delta_{D-1}$ such that $P_\sigma$ has rank $1$ for all $\sigma \in \Sigma$.
      The random variables $\{X_i: i \in \sigma\}$ are completely independent for $\sigma \in \Sigma$.
      
      This model was studied in {\em Marginal Independence Models} by Tobias Boege, Sonja Petrovi\'c and Bernd Sturmfels, 
      **INSERT CITATION TO OUR ARXIV PAPER -- CHANGE THIS LINK ** @HREF"https://arxiv.org/abs/1701.07130"@. 
      The rank 1 constraints translate to vanishing of $2\times 2$ minors of flattenings of the tensor $P$. 
      This gives the ideal of the marginal independence model in  {\em natural joint probability coordinates}. 
      Under a linear change of coordinates, the ideal becomes toric. 
      The main method, @TO marginalIndependenceIdeal@, computes this toric ideal. 
      The new ring variables are called {\em Moebius coordinates}. 
      
      In the following example, as usual, we specify a complex by its facets. 
    Example
      sigma = {{1,2},{1,3},{2,3}} 
      I = marginalIndependenceIdeal sigma
    Text
      The matrix defining this toric ideal   is given by  the following polytope: 
    Example
      P = binaryPolytope sigma 
      peek values P -- ...and  how better do i see the vertices?! 
    Text
      To appreciate the linear change  of coordinates, take a look at the map:
    Example
      phi = moebius sigma
    Text 
      The target and source rings are in natural probabilities and relevant Moebius coordinates, respectively. 
      ***SAY A FEW MORE WORDS IN DETAIL?...***
    Example
      target phi
      source phi
    Text
      This is  the 3-cycle as in Example 1 in Boege-Petrovi\'c-Sturmfels. 
      The example also appears in S. Sullivant: {\em Algebraic Statistics} (2018)  
      and G. Kirkup: {\em Random variables with completely independent subcollections} (2007). 
  SeeAlso 
    marginalIndependenceIdeal
    moebius
///
doc ///
  Key
    marginalIndependenceIdeal
    (marginalIndependenceIdeal, List)
  Headline
    toric ideal of marginal independence in Moebius coordinates 
  Usage
    marginalIndependenceIdeal(List)
  Inputs
    FF : List
      of facets of a simplicial complex S
  Outputs
    : Ideal 
      of relations on the relevant Moebius coordinates which is the prime ideal of the statistical model of 
      marginal independence encoded by S
  Description
    Text
      .... This ideal was obtained using @TO toricMarkov@ command from the 4ti2.
    Example
      FF = {{1,2},{2,3},{3,4},{4,1}}
      marginalIndependenceIdeal FF
    Text
      WE PROBABLY WANT TO DO DIFFERENT EXAMPLE HERE. 
      COMPARE TO PAPER. 
      DO A SMALL INTERESTING ONE, AND A LARGE INTERESTING ONE? 
///
doc ///
  Key
    moebius
    (marginalIndependenceIdeal, List)
  Headline
    toric ideal of marginal independence in Moebius coordinates 
  Usage
    marginalIndependenceIdeal(List)
  Inputs
    FF : List
      of facets of a simplicial complex S
  Outputs
    : RingMap 
      a linear change of coordinates from natural joint probabilities to relevant Moebius coordinates 
  Description
    Text
      .... TO BE WRITTEN. 
    Example
      FF = {{1,2},{2,3},{3,4}} 
      phi = moebius FF
      vars source phi 
      vars target phi
      --do we need an example of some variable mapping to something.. maybe later. 
///

--***************************************************************************--
--                          Tests             	                             --
--***************************************************************************--

TEST /// -- marginalIndependenceIdeal for the complex [12][23]
    sigma = {{1,2},{2,3}};
    Isigma = marginalIndependenceIdeal sigma;
    -- running some simple checks; seems like encoding teh q's by hand would not be future-version safe..
    assert (dim Isigma == 4); 
    assert (numcols mingens Isigma == 3); 
    assert (degree Isigma == 3); 
///

-- each exported method should have at least one meaningful test. 


--***************************************************************************--
end   


--***************************************************************************--
-- Useful commands
--***************************************************************************--
restart
installPackage("MarginalIndependenceModels", RemakeAllDocumentation=>true, RerunExamples=>true)
loadPackage("MarginalIndependenceModels",Reload=>true)
viewHelp "MarginalIndependenceModels"
viewHelp marginalIndependenceIdeal
installPackage("MarginalIndependenceModels")
check MarginalIndependenceModels


