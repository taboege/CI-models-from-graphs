loadPackage "GraphicalModels";
loadPackage "Matroids";

-- Generate a list of all (n choose 2)*2^(n-2) conditional independence symbols
CIstmts = method();
--------------------------------------------------------------------------------
CIstmts List := N -> (
  flatten apply(subsets(N, 2) / sort, ij -> (
    M = toList(difference(set N, ij));
    apply(subsets(M) / sort, K -> {ij#0, ij#1, K})
  ))
);

CIstmts Set := N -> CIstmts(sort(toList N));
CIstmts ZZ  := n -> CIstmts(toList(1 .. n));
--------------------------------------------------------------------------------

-- Return all CI statements satisfied by the rank function of a matroid.
-- Given a graph, take its cycle matroid.
CIrelation = method();
--------------------------------------------------------------------------------
CIrelation Matroid := M -> (
  select(CIstmts(M.groundSet), ijK -> (
    i = set({ijK#0}); j = set({ijK#1}); K = set(ijK#2);
    rank(M, K) + rank(M, i+j+K) == rank(M, i+K) + rank(M, j+K)
  ))
);

CIrelation Graph := G -> CIrelation(matroid(G));
--------------------------------------------------------------------------------

-- The Matroids package insists on using 0-based numbers for the ground set
-- of a graphic matroid. This is incompatible with GraphicalModels. Provide
-- a function for fixing this off-by-one discrepancy.
upByOne = I -> apply(I, ijK -> {
  ijK#0 + 1, ijK#1 + 1,
  apply(ijK#2, k -> k + 1)
});

-- The return value of CIrelation is not packed up enough either for
-- GraphicalModels's conditionalIndependenceIdeal. Convert {i,j,K}
-- to {{i},{j},K} (where every component is a list).
packCI = I -> apply(I, ijK -> {
  {ijK#0}, {ijK#1}, ijK#2
});

-- Return the ideal on the atomic probabilities generated by the
-- CI relation of a matroid (or graph, via its cycle matroid).
binaryIdeal = method();
--------------------------------------------------------------------------------
binaryIdeal Matroid := M -> (
  R = markovRing(#M.groundSet : 2);
  conditionalIndependenceIdeal(R, packCI(upByOne(CIrelation M)))
);

binaryIdeal Graph := G -> binaryIdeal(matroid(G));
--------------------------------------------------------------------------------

G = graph{{1,2},{2,3},{3,1}}        -- C_3
C4 = graph{{1,2},{2,3},{3,4},{4,1}}; -- C_4
--G = graph(subsets(4,2));            -- K_4
CIrelation C4/print;

C5 = graph{{1,2},{2,3},{3,4},{4,5},{5,1}};
CIrelation C5/print;
