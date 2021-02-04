EPS = 10^(-Round[MachinePrecision] - 10);
digits = 20;

$Assumptions = { p2 >= 0, m12 >= 0, m22 >= 0, m2 >= 0, q2 > 0 };

fB[x_] := Re[-1 + Log[1 - x] - x Log[1 - 1/x]]

B0[p2_, 0, 0, q2_, eps_:EPS] :=
    2 - Log[p2 / q2]

B0[0, 0, m22_, q2_, eps_:EPS] :=
    1 - Log[m22 / q2]

B0[p2_, 0, m22_, q2_, eps_:EPS] :=
    2 - Log[p2 / q2] + (m22 - p2)/p2 Log[(m22 - p2 - I eps)/m22]

B0[p2_, m12_, 0, q2_, eps_:EPS] :=
    B0[p2, 0, m12, q2, eps]

B0[0, m2_, m2_, q2_, eps_:EPS] :=
    -Log[m2 / q2]

B0[p2_, m12_, m22_, q2_, eps_:EPS] :=
    B0[p2, m22, m12, q2, eps] /; m22 > m12

B0[p2_, m12_, m22_, q2_, eps_:EPS] :=
    Module[{x, xp, xm,
            s = p2 - m22 + m12,
            imin = m12 - I eps},
        x  = Sqrt[s^2 - 4 p2 imin];
        xp = (s + Sign[s] x) / (2 p2);
        xm = imin / (xp p2);
        -Log[p2/q2] - fB[xp] - fB[xm]
    ]

(* B0[p2, 0, 0, q2] *)
B0xx[p2_, 0, q2_] := 2 - Log[p2 / q2]

(* B0[p2, 0, 0, q2] *)
B0xx[0, m2_, q2_] := -Log[m2 / q2]

(* B0[p2, m2, m2, q2] *)
B0xx[p2_, m2_, q2_] :=
    Module[{sq},
           If[p2 <= 4 m2,
              2 - Log[m2 / q2] - 2 Sqrt[4 m2 / p2 - 1] ArcSin[Sqrt[p2 / (4 m2)]],
              sq = Sqrt[1 - 4 m2 / p2];
              2 - Log[m2 / q2] + sq Log[p2 (1 - sq) / (2 m2) - 1]
           ]
    ]

(*************************************************************
   B0xx
   Format: p^2, m1^2, m2^2, q^2, B0xx[p^2, m1^2, q^2]
 *************************************************************)

MakePointB0xx[p2_, m2_, q2_] :=
    { p2, m2, m2, q2, B0xx[p2, m2, q2] }

data = ParallelMap[N[MakePointB0xx @@ #, digits]&,
    Join[
        Flatten[#,1]& @ Table[{10^n, 10^k, 1}, {n, -15, 15}, {k, -15, 15}],
        Flatten[#,3]& @ Table[{1 + s1 10^n, 1 + s2 10^k, 1},
                              {n, -15, -1}, {k, -15, -1}, {s1, -1, 1}, {s2, -1, 1}],
        {
            {1, 0, 1},
            {0, 1, 1}
        }
    ]
]

Export["B0xx.dat", data]

(*************************************************************
   B0
   Format: p^2, m1^2, m2^2, q^2, B0[p^2, m1^2, m2^2, q^2]
 *************************************************************)

MakePointB0[p2_, m12_, m22_, q2_] :=
    { p2, m12, m22, q2, B0[p2, m12, m22, q2] }

data = ParallelMap[N[MakePointB0 @@ #, digits]&,
    Join[
        Flatten[#,2]& @ Table[{10^n, 10^k, 10^l, 1}, {n, -15, 15},
                              {k, -15, 15}, {l, -15, 15}],
        Flatten[#,5]& @ Table[{1 + s1 10^n, 1 + s2 10^k, 1 + s3 10^l, 1},
                              {n, -15, -1}, {k, -15, -1}, {l, -15, -1},
                              {s1, -1, 1}, {s2, -1, 1}, {s3, -1, 1}],
        Table[{1, 10^n  , 2*10^n, 1}, {n, -15, 15}],
        Table[{1, 2*10^n, 10^n  , 1}, {n, -15, 15}],
        {
            {0, 0, 1, 1},
            {0, 1, 0, 1},
            {1, 0, 1, 1},
            {1, 1, 0, 1},
            {1, 0, 0, 1}
        }
    ]
]

Export["B0.dat", data]
