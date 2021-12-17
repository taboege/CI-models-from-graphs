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
  q1S1 - z*a1*(1-b1-b2)*c1,
  q1S2 - z*a1*(1-b1-b2)*c2,
  q2S1 - z*a2*(1-b1-b2)*c1,
  q2S2 - z*a2*(1-b1-b2)*c2,
  qS11 - z*(1-a1-a2)*b1*c1,
  qS12 - z*(1-a1-a2)*b1*c2,
  qS21 - z*(1-a1-a2)*b2*c1,
  qS22 - z*(1-a1-a2)*b2*c2,
  q1SS - z*a1*(1-b1-b2)*(1-c1-c2),
  q2SS - z*a2*(1-b1-b2)*(1-c1-c2),
  qS1S - z*(1-a1-a2)*b1*(1-c1-c2),
  qS2S - z*(1-a1-a2)*b2*(1-c1-c2),
  qSS1 - z*(1-a1-a2)*(1-b1-b2)*c1,
  qSS2 - z*(1-a1-a2)*(1-b1-b2)*c2,
  qSSS - z*(1-a1-a2)*(1-b1-b2)*(1-c1-c2)
);

-- Implicitize
M = eliminate({z, a1, b1, c1, a2, b2, c2}, I);

<< betti mingens M << endl << endl;
<< betti res M << endl << endl;
exit();
