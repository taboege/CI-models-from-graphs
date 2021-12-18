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
    "moebiusRing",
    "Coefficients",
    "VariableName",
    "RelevantMoebius",
    -- "quadraticModel", 
    "marginalIndependenceIdeal"
    }

needsPackage "FourTiTwo"; --we compute a Markov basis
needsPackage "Polyhedra"; --we construct a polytope

-- the IndexedVariable names for two coordinate rings: 
p = local p;
q = local q; -- this will go away after moebiusRing is finalized. 

moebiusRingData = local moebiusRingData
moebiusVariables = local moebiusVariables


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


------------------------------------------------------------------------------------------------
-- moebiusRing List
-* 
*****this is (almost) literarlly doing what GraphicalModels.m2 does for MarkovRing *****
*- 
-- Outputs a polynomial ring whose indeterminates are Moebius coordinates for discrete
-- random variables with a given number of states. 
-- d should be a sequence of integers di >= 1
--
-- NOTE: there is a mutable hash table of all Moebius rings created, so as to not re-create rings!
-- the hashtable is indexed by the sequence d, the coefficient ring kk, and the variable name p, 
-- as this information identifies the Moebius ring uniquely. 
------------------------------------------------------------------------------------------------

toSymbol = (p) -> (
     if instance(p,Symbol) then p
     else
     if instance(p,String) then getSymbol p
     else
     error ("expected a string or symbol, but got: ", toString p))


moebiusRingList := new MutableHashTable;

-- INPUT
-- FF, simplicial complex (as a list of facets) 
-- FUTURE INPUT
--  one day we may want to do a general state ring, so in that case, 
--  d = sequence of integers (d_1,...,d_n) for the number of states of the variables 1,...,n.
--  today, we only use binary models. 
--  so we preset d = (2,...,2). 
-- OUTPUT
--  poly ring whose variables are all moebius labels. 
-- optional input: 
-- RelevantMoebius  = whetherto include ONLY relevant mob.coords. useful wehn you just want ideal and not whole map...
moebiusRing = method(Dispatch=>Thing, Options=>{Coefficients=>QQ, VariableName=>"q", RelevantMoebius=>true})
--moebiusRing Sequence := Ring => opts -> (FF,d) -> (
moebiusRing List := Ring => opts -> FF -> (
    --if any(d, di -> not instance(di,ZZ) or di <= 0)
    --    then error "moebiusRing expected nonnegative integers";
     d := (max flatten FF : 2);  --- don't forget that indexing starts at 0. 
     kk := opts.Coefficients;
     q := toSymbol opts.VariableName;
     if not moebiusRingList#?(FF,d,kk,q) then (
	 --start := (#d):0;
	 --vlist := start .. d-1;
	 N := sort unique flatten FF;  
	 FF = complexify FF; 
	 -- if opts.RelevantMoebius then R := kk(monoid [ apply(relevantMoebiusLabels(FF), A -> q_A) ]) else R := kk(monoid [ apply(moebiusLabels(N), A -> q_A) ]);  
	 if opts.RelevantMoebius then (
	  --R := kk(monoid [ select(subsets(N), A -> member(A, FF)) / sort / (A -> q_A)]); 
	  R := kk(monoid [ apply(relevantMoebiusLabels(FF), A -> q_A) ]); 
	  H := new HashTable from apply(relevantMoebiusLabels(N), A -> A=>q_A); 
  	  R.moebiusVariables = H; 
	 ) else (
	  R = kk(monoid [ apply(moebiusLabels(N), A -> q_A) ]);
	  H = new HashTable from apply(moebiusLabels(N), A -> A=>q_A); 
  	  R.moebiusVariables = H; 
	 );

         R.moebiusRingData = d;
	 moebiusRingList#(FF,d,kk,q) = R;
	 );	  
    moebiusRingList#(FF,d,kk,q)
    ); 


----------------------------------------
--      Moebius setup
----------------------------------------
-*
Input:
  FF = simplicial complex, given as a list of facets
Output:
  The linear ring map which transforms the natural probability variables p into Moebius coordinates q. 
*-
moebius = method(
    TypicalValue => List
    );
moebius List := (PolynomialRing,PolynomialRing,RingMap) => FF -> (
  q = symbol q;
  p = symbol p;
  N := sort unique flatten FF;
  Q := QQ[ apply(moebiusLabels(N), A -> q_A) ];
  
  --  Q := moebiusRing( FF, RelevantMoebius=>false); -- we need the entire ring here, not just relevant q's!! 
  P := QQ[ apply(binaryLabels(N),  a -> p_a) ];  
  
  mu := map(P, Q, apply(moebiusLabels(N), A -> q_A => (
  --mu := map(P, Q, apply(Q.moebiusLabels, A -> q_A => (
  --mu := map(P, Q, apply(keys Q.moebiusVariables, A -> (Q.moebiusVariables)#A => (
    -- D is the set where 0/1 should vary in the summation
    D := select(N, d -> not member(d, A));
    sum apply(apply(subsets(D), B -> (
      a := new MutableList from #N:0;
      apply(D, d -> (
        a#(position(N, i -> i == d)) = if member(d, B) then 1 else 0;
      ));
      toSequence a
    )), a -> p_a)
  )));

  (P,Q,mu)
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
--  q = symbol q; -- this is guided by how indexed variables were handled in GraphicalModels.m2 :) 
  N := sort unique flatten FF;
  B := binaryPolytope FF;
  FF = complexify FF;
  Q := moebiusRing(FF,RelevantMoebius =>true);-- QQ[ select(subsets(N), A -> member(A, FF)) / sort / (A -> q_A) ];
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
-*
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
*-

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
moebiusLabels List := List => N -> (
    subsets(N)
    ); 
-*
moebiusLabels = method(Dispatch=>Thing, Options=>{RelevantMoebius=>true})
moebiusLabels List := List => opts -> N -> (
    if opts.RelevantMoebius then ( 
	-- THIS ONE NEEDS FF AND N;
	select(subsets(N), A -> member(A, FF)) / sort / (A -> q_A)
	) else subsets(N) -- THIS ONE ONLY N.
);
*-
relevantMoebiusLabels = method( 
    TypicalValue => List
    );
relevantMoebiusLabels List := List => FF -> (
    N := sort unique flatten FF;
    FF = complexify FF;
    select(subsets(N), A -> member(A, FF)) / sort / (A -> A)
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
      
      This model was studied in {\em Marginal Independence Models} by Tobias Boege, Sonja Petrovi\'c and Bernd Sturmfels, @HREF"https://arxiv.org"@. 
      The rank 1 constraints translate to vanishing of $2\times 2$ minors of flattenings of the tensor $P$. 
      This gives the ideal of the marginal independence model in  {\em natural joint probability coordinates}. 
      Under a linear change of coordinates, the ideal becomes toric. 
      The main method, @TO marginalIndependenceIdeal@, computes this toric ideal. 
      The new ring variables are called {\em Moebius coordinates}. 
      
      In the following example, as usual, we specify a complex by its facets. 
    Example
      sigma = {{1,2},{1,3},{2,3}} 
      moebiusRing sigma
      I = marginalIndependenceIdeal sigma
    Text
      The matrix defining this toric ideal   is given by  the following polytope: 
    Example
      needsPackage "Polyhedra"; 
      P = binaryPolytope sigma 
      vertices P
    Text
      To appreciate the linear change  of coordinates, take a look at the map:
    Example
      (P,Q,mu) = moebius sigma
    Text 
      The target and source rings are in natural probabilities and relevant Moebius coordinates, respectively. 
      ***SAY A FEW MORE WORDS IN DETAIL?...***
    Example
      target mu
      source mu
      vars P
      vars Q
    Text
      This is  the 3-cycle as in Example 1 in Boege-Petrovi\'c-Sturmfels. 
      The example also appears in S. Sullivant: {\em Algebraic Statistics} (2018)  
      and G. Kirkup: {\em Random variables with completely independent subcollections} (2007).
  Caveat
    It  is worth noting that in this package, the  ambient polynomial ring uses only  the relevant Moebius coordinates,
    while the full model lives in a larger space. Hence there will be  a discrepancy 
    comparing the theoretical model results and properties of the ideal $I$, for example, 
    dimension is offset by the number of irrelevant Moebius coordinates. 
    The true (affine) dimension of the model is $2^n$ - codim I. 
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
    binaryPolytope
    (binaryPolytope, List)
  Headline
    the polytope of the marginal independence model for binary random variables
  Usage
    binaryPolytope(List)
  Inputs
    FF : List
      of facets of a simplicial complex S
  Outputs
    : Polyhedron 
      with vertices 
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
    (moebius, List)
  Headline
    toric ideal of marginal independence in Moebius coordinates 
  Usage
    moebius(List)
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
      (P,Q,mu) = moebius FF
      vars source mu 
      transpose vars target mu
      (inverse mu) (vars P)_0
      mu (vars Q)_0
      --do we need an example of some variable mapping to something.. yes. a binomial mapping to a hot mess :) :)  
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


