needsPackage "FourTiTwo";
needsPackage "Polyhedra";
needsPackage "PhylogeneticTrees";

-- Generate the simplicial complex from a sequence of sets
complexify = FF -> unique flatten apply(FF, F -> subsets(F) / sort);

-- Conditional independence statements corresponding to the
-- independences in the complex FF.
CIstmts = N -> (
  flatten apply(subsets(N, 2) / sort, ij -> (
    M := toList(difference(set N, ij));
    apply(subsets(M) / sort, K -> {ij#0, ij#1, K})
  ))
);
indepCI = FF -> unique flatten apply(FF, CIstmts);

-- The 0/1 polytope whose vertices are the faces of FF.
binaryPolytope = FF -> (
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

binaryLabels = N -> (
  A := set{0,1};
  toList(fold((X,Y) -> (X ** Y) / splice, #N:A)) / toSequence
);
moebiusLabels = N -> (
  subsets(N)
);

-- Get the whole Moebius setup: the ring in p-variables,
-- the ring in q-variables and the linear map.
moebius = FF -> (
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

-- The Segre map in P coordinates.
segre = (FF, P) -> (
  N := sort unique flatten FF;
  S := QQ[ flatten apply(N, i -> { p_(i,0), p_(i,1) }) ];
  map(S, P, apply(binaryLabels(N), a -> p_a => (
    product apply(N, i -> (
      v := a#(position(N, j -> j == i));
      p_(i,v)
    ))
  )))
);

-- The linear space of margin-zero tensors in P coordinates.
marg = (FF, P) -> (
  N := sort unique flatten FF;
  ideal(flatten apply(FF / sort, F -> (
    D := select(N, d -> not member(d, F));
    apply(subsets(F), G -> (
      sum apply(apply(subsets(D), B -> (
        a := new MutableList from #N:0;
        BG := set(B) + set(G);
        apply(N, i -> (
          a#(position(N, j -> j == i)) = if member(i, BG) then 1 else 0;
        ));
        toSequence a
      )), a -> p_a)
    ))
  )))
);

-- The join of marg and segre for a given simplicial complex.
binaryJoinModel = FF -> (
  f := moebius FF;
  P := target f;
  S := kernel segre(FF, P);
  L := marg(FF, P);
  joinIdeal(S, L)
);

-- The model according to Theorem 5 in [Nov 28]. 
-- the special quadrics only -  we know this isn't the full ideal,
-- but maybe these are conjectured to generate ideal for matroids (White).
binaryModel = FF -> (
  N := sort unique flatten FF;
  f := moebius(FF);
  Q := source f; use Q;
  FF = complexify(FF) / set;
  good = (s,t) -> member(set(s) + set(t), FF);
  ideal flatten apply(subsets(N), s -> (
    select(unique apply(select(subsets(N), t -> good(s,t)), t -> (
      i = sort toList(set(s) * set(t));
      u = sort toList(set(s) + set(t));
      q_s * q_t - q_i * q_u
    )), f -> f != 0)
  ))
);

polytopeModel = FF -> (
  N := sort unique flatten FF;
  B := binaryPolytope FF;
  FF = complexify FF;
  Q := QQ[ select(subsets(N), A -> member(A, FF)) / sort / (A -> q_A) ];
  --toricGroebner(Polyhedra$vertices(B), Q)
  toricMarkov(Polyhedra$vertices(B), Q)
);



