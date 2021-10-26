-- Three binary random variables with facets {1,2}, {1,3}, {2,3},
-- as in Kirkup's paper. The 0 index denotes "bullet" or "summation".
S = QQ[
  x110, x120, x210, x220,
  x101, x102, x201, x202,
  x011, x012, x021, x022,
  x100, x200,
  x010, x020,
  x001, x002,
  x000
];

-- The ring with all atomic probabilities, codomain of the map
-- defining K.
R = QQ[p111, p112, p121, p122, p211, p212, p221, p222];

-- The ring corresponding to factorizations, codomain of the map
-- defining Q.
T = QQ[y10, y11, y12, y20, y21, y22, y30, y31, y32];

-- Summation along the zeros.
K = kernel map(R, S, {
  p111 + p112, p121 + p122, p211 + p212, p221 + p222,
  p111 + p121, p112 + p122, p211 + p221, p212 + p222,
  p111 + p211, p112 + p212, p121 + p221, p122 + p222,
  p111 + p112 + p121 + p122, p211 + p212 + p221 + p222,
  p111 + p112 + p211 + p212, p121 + p122 + p221 + p222,
  p111 + p121 + p211 + p221, p112 + p122 + p212 + p222,
  p111 + p112 + p121 + p122 + p211 + p212 + p221 + p222
});

-- Products along the ground set.
Q = kernel map(T, S, {
  y11*y21*y30, y11*y22*y30, y12*y21*y30, y12*y22*y30,
  y11*y20*y31, y11*y20*y32, y12*y20*y31, y12*y20*y32,
  y10*y21*y31, y10*y21*y32, y10*y22*y31, y10*y22*y32,
  y11*y20*y30, y12*y20*y30,
  y10*y21*y30, y10*y22*y30,
  y10*y20*y31, y10*y20*y32,
  y10*y20*y30
});

print dim(S/K) -- expected 7
print dim(S/Q) -- unexpected 7 (should be 6 according to the paper)
