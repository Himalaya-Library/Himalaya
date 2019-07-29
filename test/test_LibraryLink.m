passedTests = 0;
failedTests = 0;

TestEqual[a_, b_] :=
    If[a === b,
       Print["Test passed: ", a, " === ", b];
       passedTests++
       ,
       Print["Test failed: ", a, " =!= ", b];
       failedTests++
      ];

TestClose[a_List, b_List, eps_:10^-15] :=
    TestClose[Sequence @@ #, eps]& /@ MapThread[List, {a, b}];

TestClose[a_, b_, eps_:10^-15] :=
    If[Abs[(a - b)] < eps || (Abs[a] > eps && Abs[(a - b)/a] < eps),
       Print["Test passed: ", a, " ~ ", b, " within ", eps];
       passedTests++
       ,
       Print["Test failed: ", a, " !~ ", b, " within ", eps];
       failedTests++
       ,
       Print["Test failed: Could not evaluate expression to bool: ", Abs[a - b] < eps];
       failedTests++
      ];

TestPoint[point_, name_] :=
    Module[{output, eps},
           output = HimalayaCalculateDMh3L[Sequence @@ (input /. point)];
           eps = precision /. point;

           Print[];
           Print["==== Testing Point ", name, " ================================"];
           Print[];

           TestEqual[hierarchyID                                    /. output, hierarchyID                                    /. (expectedOutput /. point)];
           TestEqual[hierarchyName                                  /. output, hierarchyName                                  /. (expectedOutput /. point)];
           TestClose[MstopMDRPrime                                  /. output, MstopMDRPrime                                  /. (expectedOutput /. point), eps];
           TestClose[MsbottomMDRPrime                               /. output, MsbottomMDRPrime                               /. (expectedOutput /. point), eps];
           TestClose[Mh2                                            /. output, Mh2                                            /. (expectedOutput /. point), eps];
           TestClose[Mh2ShiftDRbarPrimeToMDRPrime                   /. output, Mh2ShiftDRbarPrimeToMDRPrime                   /. (expectedOutput /. point), eps];
           TestClose[Mh2ShiftDRbarPrimeToH3m                        /. output, Mh2ShiftDRbarPrimeToH3m                        /. (expectedOutput /. point), eps];
           TestClose[expansionUncertainty                           /. output, expansionUncertainty                           /. (expectedOutput /. point), eps];
           TestClose[Mh2EFT                                         /. output, Mh2EFT                                         /. (expectedOutput /. point), eps];
           TestClose[Mh2FO                                          /. output, Mh2FO                                          /. (expectedOutput /. point), eps];
           TestClose[Mh2FOAt                                        /. output, Mh2FOAt                                        /. (expectedOutput /. point), eps];
           TestClose[lambda                                         /. output, lambda                                         /. (expectedOutput /. point), eps];
           TestClose[lambdaUncertainty                              /. output, lambdaUncertainty                              /. (expectedOutput /. point), eps];
           TestClose[lambdaShiftDRbarPrimeToMSbar                   /. output, lambdaShiftDRbarPrimeToMSbar                   /. (expectedOutput /. point), eps];
          ];

(* load LibrayLink *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];

InitializeHimalaya[libPath];

HimalayaInfoMessage[s_]  := Print[s];
HimalayaErrorMessage[s_] := Print[s];

Off[General::stop];

(*
himalaya::Parameters test_point() {
   himalaya::Parameters pars;

   const double MS = 2000.;
   const double MS2 = MS*MS;
   const double Xt = sqrt(6.)*MS;
   const double tb = 20.;

   pars.scale = MS;
   pars.mu = MS;
   pars.g1 = 0.46;
   pars.g2 = 0.65;
   pars.g3 = 1.166;
   pars.vd = 246*std::cos(std::atan(tb));
   pars.vu = 246*std::sin(std::atan(tb));
   pars.mq2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.md2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.mu2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.ml2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.me2 << MS2, 0, 0,
               0, MS2, 0,
               0, 0, MS2;
   pars.Au(2,2) = Xt + pars.mu/tb;
   pars.Yu(2,2) = 0.862;
   pars.Yd(2,2) = 0.133;
   pars.Ye(2,2) = 0.101;
   pars.MA = MS;
   pars.M1 = MS;
   pars.M2 = MS;
   pars.MG = MS;

   return pars;
}
*)

MakePoint[MS_?NumericQ, TB_?NumericQ, Xt_?NumericQ] := {
    input -> {
        settings -> {
            bottom -> False,
            verbose -> True
        },
        parameters -> {
            scale -> MS,
            mu -> MS,
            g1 -> 0.46,
            g2 -> 0.65,
            g3 -> 1.166,
            vd -> 246*Cos[ArcTan[TB]],
            vu -> 246*Sin[ArcTan[TB]],
            mq2 -> MS^2 IdentityMatrix[3],
            md2 -> MS^2 IdentityMatrix[3],
            mu2 -> MS^2 IdentityMatrix[3],
            ml2 -> MS^2 IdentityMatrix[3],
            me2 -> MS^2 IdentityMatrix[3],
            Au -> {{0,0,0},
                   {0,0,0},
                   {0,0, Xt + MS/TB }},
            Ad -> 0 IdentityMatrix[3],
            Ae -> 0 IdentityMatrix[3],
            Yu -> {{0,0,0},
                   {0,0,0},
                   {0,0, 0.862 }},
            Yd -> {{0,0,0},
                   {0,0,0},
                   {0,0, 0.133 }},
            Ye -> {{0,0,0},
                   {0,0,0},
                   {0,0, 0.101 }},
            MA -> MS,
            M1 -> MS,
            M2 -> MS,
            M3 -> MS
        }
    },
    expectedOutput -> {
        hierarchyID -> 1, hierarchyName -> h32q2g, 
        MstopMDRPrime -> {1807.421718155926, 2176.2140733830706}, 
        Mh2 -> {
            {{3.9900456677813977*^6, -199915.84939351628}, {-199915.84939351628, 18267.112558603498}},
            {{-639.5967086048405, 38.1108098347363}, {38.1108098347363, 10354.694221748936}},
            {{-2.067756466175066, 47.44906816768487}, {47.44906816768487, 1872.6875317790061}},
            {{-4.1862855658588956, 26.44030520827635}, {26.44030520827635, 695.9601307370843}}
        },
        Mh2ShiftDRbarPrimeToMDRPrime -> {
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{-0.06556209414766911, 6.451828360749449}, {6.451828360749449, -47.47567025934484}}
        },
        Mh2ShiftDRbarPrimeToH3m -> {
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{-1.5317742384518245, 1.9557280682757228}, {1.9557280682757228, 7.446661229467915}}
        },
        expansionUncertainty -> {0., 0., 0.29936990221850124, 0.029834152694507742}, 
        Mh2EFT -> {8230.066622149989, 8168.571812031212, 830.6875834664329, 695.6106523929473},
        Mh2FO -> {8229.896089958725, 8148.2372899382335, 819.1872548102838, 696.8358059760103},
        Mh2FOAt -> {0., 44653.01081718066, 820.6212844284601, 696.8358059760103},
        lambda -> {0.13599819257964818, 0.06251358668613667, 0.0014909887921878413, 0.00031561329131293657}, 
        lambdaUncertainty -> {0., 0., 0., 0.002031176262211924},
        lambdaShiftDRbarPrimeToMSbar -> {0., 0., 7.2898259577239126*^-6, -0.000771028340668861}
    },
    precision -> 1*^-5
};

TestPoint[MakePoint[2000, 20, Sqrt[6]*2000], "MS=2000, TB=20, Xt=Sqrt[6]*MS"];

Print["Number of passed tests: ", passedTests];
Print["Number of failed tests: ", failedTests];

Quit[failedTests];
