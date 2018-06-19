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
           TestClose[lambda                                         /. output, lambda                                         /. (expectedOutput /. point), eps];
           TestClose[lambdaUncertainty                              /. output, lambdaUncertainty                              /. (expectedOutput /. point), eps];
           TestClose[lambdaShiftDRbarPrimeToMSbar                   /. output, lambdaShiftDRbarPrimeToMSbar                   /. (expectedOutput /. point), eps];
           TestClose[deltaLambda3LoopH3m                            /. output, deltaLambda3LoopH3m                            /. (expectedOutput /. point), eps];
           TestClose[deltaLambda3LoopH3mShiftDRbarPrimeToMSbar      /. output, deltaLambda3LoopH3mShiftDRbarPrimeToMSbar      /. (expectedOutput /. point), eps];
           TestClose[deltaLambda3LoopEFT                            /. output, deltaLambda3LoopEFT                            /. (expectedOutput /. point), eps];
           TestClose[deltaLambda3LoopEFTShiftDRbarPrimeToMSbar      /. output, deltaLambda3LoopEFTShiftDRbarPrimeToMSbar      /. (expectedOutput /. point), eps];
          ];

(* load LibrayLink *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];

InitializeHimalaya[libPath];

HimalayaInfoMessage[s_]  := Print[s];
HimalayaErrorMessage[s_] := Print[s];

Off[General::stop];

(* SPS2 *)
pointSPS2 = {
    input -> {
        settings -> {
            bottom -> False,
            verbose -> True
        },
        parameters -> {
            scale -> 1.11090135*^03,
            mu -> 3.73337018*^02,
            g3 -> 1.06187116*^00,
            vd -> 2.51008404*^01,
            vu -> 2.41869332*^02,
            mq2 -> { {2.36646981*^06, 0, 0},
                     {0, 2.36644973*^06, 0},
                     {0, 0, 1.63230152*^06} },
            md2 -> { {2.35612778*^06, 0, 0},
                     {0, 2.35610884*^06, 0},
                     {0, 0, 2.31917415*^06} },
            mu2 -> { {2.35685097*^06, 0, 0},
                     {0, 2.35682945*^06, 0},
                     {0, 0, 9.05923409*^05} },
            Ad -> { {0, 0, 0},
                    {0, 0, 0},
                    {0, 0, -784.3356416708631} },
            Au -> { {0, 0, 0},
                    {0, 0, 0},
                    {0, 0, -527.8746242245387} },
            MA -> 1.48446235*^03,
            M3 -> 6.69045022*^02,
            MW -> 8.04001915*^01,
            MZ -> 8.97608307*^01,
            Mt -> 1.47685846*^02,
            Mb -> 2.38918959*^00,
            MSt -> { 9.57566721*^02, 1.28878643*^03 },
            MSb -> { 1.27884964*^03, 1.52314587*^03 },
            s2t -> Sin[2*ArcSin[1.13197339*^-01]],
            s2b -> Sin[2*ArcSin[-9.99883015*^-01]]
        }
    },
    expectedOutput -> {
        hierarchyID                                    -> 2,
        hierarchyName                                  -> h3q22g,
        MstopMDRPrime                                  -> {955.6409734756253, 1287.5049731244299},
        Mh2                                            -> {
            {{2.18023*^06, -227080}, {-227080, 31451.3}},
            {{-2.92153, 100.252}, {100.252, 5277.06}},
            {{-0.307806, 17.7881}, {17.7881, 1256.76}},
            {{-0.0143744, 5.62494}, {5.62494, 300.795}}
        },
        Mh2ShiftDRbarPrimeToMDRPrime                   -> {
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{0.00487318, 0.0179684}, {0.0179684, 0.738644}}
        },
        Mh2ShiftDRbarPrimeToH3m                        -> {
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{0, 0}, {0, 0}},
            {{-0.0037627, 0.0622916}, {0.0622916, -0.631601}}
        },
        expansionUncertainty                           -> {0, 0, 0.0378711, 0.0243905},
        Mh2EFT                                         -> {7717.27, 5241.38, 1245.02, 303.339},
        lambda                                         -> {0.13051177168737485, 0.00518866, 0.000211537, -1.17351*^-05},
        lambdaUncertainty                              -> {0, 0, 0, 0.000102286},
        lambdaShiftDRbarPrimeToMSbar                   -> {0, 0, -3.9067*^-07, -4.62254*^-05},
        deltaLambda3LoopH3m                            -> {8.33505*^-06, 0.00012236},
        deltaLambda3LoopH3mShiftDRbarPrimeToMSbar      -> -4.62215*^-05,
        deltaLambda3LoopEFT                            -> {-1.17351*^-05, 0.000102286},
        deltaLambda3LoopEFTShiftDRbarPrimeToMSbar      -> -4.62254*^-05
    },
    precision -> 1*^-5
};

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
            {{-4.186285565859146, 26.440305208276936}, {26.440305208276936, 695.9601307370843}}
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
        Mh2EFT -> {8230.066622149989, 10331.078614427297, 1881.1471720931024, 689.1169793778004},
        lambda -> {0.13599819257964818, 0.06251358668613667, 0.0014909887921878413, 0.00021462203207113314}, 
        lambdaUncertainty -> {0., 0., 0., 0.0010928804743410628},
        lambdaShiftDRbarPrimeToMSbar -> {0., 0., 7.2898259577239126*^-6, -0.000771028340668861}, 
        deltaLambda3LoopH3m -> {0.000226477371058642, 0.0012078034614150654}, 
        deltaLambda3LoopH3mShiftDRbarPrimeToMSbar -> -0.0006679606925823672, 
        deltaLambda3LoopEFT -> {0.00021462203207113314, 0.0010928804743410628}, 
        deltaLambda3LoopEFTShiftDRbarPrimeToMSbar -> -0.000771028340668861
    },
    precision -> 1*^-5
};

TestPoint[pointSPS2, "SPS2"];
TestPoint[MakePoint[2000, 20, Sqrt[6]*2000], "MS=2000, TB=20, Xt=Sqrt[6]*MS"];

Print["Number of passed tests: ", passedTests];
Print["Number of failed tests: ", failedTests];

Quit[failedTests];
