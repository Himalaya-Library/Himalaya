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

           TestEqual[renormalizationScheme                          /. output, renormalizationScheme                          /. (expectedOutput /. point)];
           TestEqual[hierarchyID                                    /. output, hierarchyID                                    /. (expectedOutput /. point)];
           TestEqual[hierarchyName                                  /. output, hierarchyName                                  /. (expectedOutput /. point)];
           TestClose[Mstop                                          /. output, Mstop                                          /. (expectedOutput /. point), eps];
           TestClose[Mh2Tree                                        /. output, Mh2Tree                                        /. (expectedOutput /. point), eps];
           TestClose[Mh21Loop                                       /. output, Mh21Loop                                       /. (expectedOutput /. point), eps];
           TestClose[Mh22Loop                                       /. output, Mh22Loop                                       /. (expectedOutput /. point), eps];
           TestClose[Mh23Loop                                       /. output, Mh23Loop                                       /. (expectedOutput /. point), eps];
           TestClose[Mh23LoopShiftDRbarPrimeToMDRPrime              /. output, Mh23LoopShiftDRbarPrimeToMDRPrime              /. (expectedOutput /. point), eps];
           TestClose[Mh23LoopShiftDRbarPrimeToH3m                   /. output, Mh23LoopShiftDRbarPrimeToH3m                   /. (expectedOutput /. point), eps];
           TestClose[expansionUncertainty                           /. output, expansionUncertainty                           /. (expectedOutput /. point), eps];
           TestClose[deltaLambda3LoopH3mDRbarPrime                  /. output, deltaLambda3LoopH3mDRbarPrime                  /. (expectedOutput /. point), eps];
           TestClose[deltaLambda3LoopH3mShiftDRbarPrimeToMSbar      /. output, deltaLambda3LoopH3mShiftDRbarPrimeToMSbar      /. (expectedOutput /. point), eps];
           TestClose[deltaLambda3LoopEFTDRbarPrime                  /. output, deltaLambda3LoopEFTDRbarPrime                  /. (expectedOutput /. point), eps];
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
        renormalizationScheme                          -> DRbarPrime,
        hierarchyID                                    -> 2,
        hierarchyName                                  -> h3q22g,
        Mstop                                          -> {957.567, 1288.79},
        Mh2Tree                                        -> {{2.18023*^06, -227080}, {-227080, 31451.3}},
        Mh21Loop                                       -> {{-2.92153, 100.252}, {100.252, 5277.06}},
        Mh22Loop                                       -> {{-0.307806, 17.7881}, {17.7881, 1256.76}},
        Mh23Loop                                       -> {{-0.0143744, 5.62494}, {5.62494, 300.795}},
        Mh23LoopShiftDRbarPrimeToMDRPrime              -> {{0.00487318, 0.0179684}, {0.0179684, 0.738644}},
        Mh23LoopShiftDRbarPrimeToH3m                   -> {{-0.0037627, 0.0622916}, {0.0622916, -0.631601}},
        expansionUncertainty                           -> {0, 0, 0.0378711, 0.0243905},
        deltaLambda3LoopH3mDRbarPrime                  -> {8.33505*^-06, 0.00012236},
        deltaLambda3LoopH3mShiftDRbarPrimeToMSbar      -> -4.62215*^-05,
        deltaLambda3LoopEFTDRbarPrime                  -> {-1.17351*^-05, 0.000102286},
        deltaLambda3LoopEFTShiftDRbarPrimeToMSbar      -> -4.62254*^-05
    },
    precision -> 1*^-5
};

TestPoint[pointSPS2, "SPS2"];

Print["Number of passed tests: ", passedTests];
Print["Number of failed tests: ", failedTests];

Quit[failedTests];
