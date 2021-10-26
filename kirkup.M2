-- Generate the simplicial complex from a sequence of sets
complexify = FF -> unique flatten apply(FF, F -> subsets(F));

-- Conditional independence statements corresponding to the
-- independences in the complex FF.
CIstmts = N -> (
  flatten apply(subsets(N, 2) / sort, ij -> (
    M = toList(difference(set N, ij));
    apply(subsets(M) / sort, K -> {ij#0, ij#1, K})
  ))
);
indepCI = FF -> unique flatten apply(FF, CIstmts);

-- Compatibility with GraphicalModels.
packCI = I -> apply(I, ijK -> {
  {ijK#0}, {ijK#1}, ijK#2
});

-- The polynomial ring in all atomic probabilities and the
-- Delta-independence ideal.
loadPackage "GraphicalModels";
binaryRing = N -> markovRing(#N : 2);
indepIdeal = FF -> (
  N = sort unique flatten FF;
  conditionalIndependenceIdeal(binaryRing(N), packCI indepCI(FF))
);

-- The ring S_\Delta in Kirkup's paper for binary random variables.
-- We use 0 to denote summation (where Kirkup uses a bullet in his
-- paper). The binary states are 1 and 2 (instead of 0 and 1),
-- just as in markovRing from GraphicalModels.
kirkLabels = FF -> (
  N = sort unique flatten FF;
  FF = complexify FF;
  flatten apply(FF, F -> (
    -- Instantiate one label for each binary vector over F.
    apply(subsets(F), A -> (
      -- The binary vector corresponding to A.
      B = apply(F, i -> if member(i, A) then 2 else 1);
      -- Splice it into the positions indexed by F, while
      -- the other positions are 0.
      a = new MutableList from #N : 0;
      apply(F, f -> (
        a#(position(N, i -> i == f)) = B#(position(F, i -> i == f));
      ));
      toSequence a
    ))
  ))
);
kirkRing = FF -> (
  vv = apply(kirkLabels(FF), a -> x_a);
  QQ[vv]
);

-- The ideal K_\Delta.
kirkMarg = FF -> (
  N = sort unique flatten FF;
  S = kirkRing FF;
  T = binaryRing N;
  kernel map(T, S, apply(kirkLabels(FF), a -> (
    -- Get all binary vectors indexed by D, the complement of F.
    D = select(N, d -> a#(position(N, i -> i == d)) == 0);
    sum apply(apply(subsets(D), A -> (
      B = apply(D, i -> if member(i, A) then 2 else 1);
      -- Splice the vector b into the zeros of a to obtain a
      -- variable from T.
      ab = new MutableList from a;
      apply(D, d -> (
        ab#(position(N, i -> i == d)) = B#(position(D, i -> i == d));
      ));
      toSequence ab
    )), a -> p_a)
  )))
);

-- The codomain of the map defining Q_\Delta.
margLabels = N -> (
  flatten apply(N, i -> (
    apply({0,1,2}, j -> toSequence {i, j})
  ))
);
margRing = N -> (
  vv = apply(margLabels(N), a -> y_a);
  QQ[vv]
);

-- The ideal Q_\Delta.
kirkIndep = FF -> (
  N = sort unique flatten FF;
  S = kirkRing FF;
  T = margRing N;
  kernel map(T, S, apply(kirkLabels(FF), a -> (
    product apply(apply(N, i -> (
      toSequence {i, a#(position(N, j -> i == j))}
    )), ij -> y_ij)
  )))
);

-- The ideal P_\Delta = radical(K_\Delta + Q_\Delta),
-- conjectured not to need the radical there.
kirkPrime = FF -> (
  S = kirkRing FF;
  radical(sub(kirkMarg(FF), S) + sub(kirkIndep(FF), S))
);
