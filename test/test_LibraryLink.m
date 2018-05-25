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

(* load LibrayLink *)
(* TODO: use better path *)
Get[FileNameJoin[{"..", "..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];
InitializeHimalaya[libPath];

Off[General::stop];

output = HimalayaCalculateDMh3L[
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
        Ab -> -784.3356416708631,
        At -> -527.8746242245387,
        MA -> 1.48446235*^03,
        MG -> 6.69045022*^02,
        MW -> 8.04001915*^01,
        MZ -> 8.97608307*^01,
        Mt -> 1.47685846*^02,
        Mb -> 2.38918959*^00,
        MSt -> { 9.57566721*^02, 1.28878643*^03 },
        MSb -> { 1.27884964*^03, 1.52314587*^03 },
        s2t -> Sin[2*ArcSin[1.13197339*^-01]],
        s2b -> Sin[2*ArcSin[-9.99883015*^-01]]
    }
];

eps = 10*^-6;

TestEqual[renormalizationScheme /. output, DRbarPrime, eps];
TestEqual[hierarchyID           /. output, 2, eps];
TestEqual[hierarchyName         /. output, h3q22g, eps];
TestClose[Mstop                 /. output, {957.567, 1288.79}, eps];
TestClose[Mh2Tree               /. output, {{2.18023*^06, -227080}, {-227080, 31451.3}}, eps];
TestClose[Mh21Loop              /. output, {{-2.92153, 100.252}, {100.252, 5277.06}}, eps];
TestClose[Mh22Loop              /. output, {{-0.307806, 17.7881}, {17.7881, 1256.76}}, eps];
TestClose[Mh23Loop              /. output, {{-0.0213027, 5.63268}, {5.63268, 299.769}}, eps];
TestClose[expansionUncertainty  /. output, {0, 0, 0.0378711, 0.0121494}, eps];
TestClose[deltaLambda3LoopHimalayaDRbarPrime             /. output, -2.3901*^-05, eps];
TestClose[deltaLambda3LoopHimalayaShiftDRbarPrimeToMSbar /. output, -4.62215*^-05, eps];
TestClose[deltaLambda3LoopHimalayaUncertainty            /. output, 1.85817*^-05, eps];
TestClose[deltaLambda3LoopEFTDRbarPrime                  /. output, -4.04954*^-05, eps];
TestClose[deltaLambda3LoopEFTShiftDRbarPrimeToMSbar      /. output, -4.62254*^-05, eps];
TestClose[deltaLambda3LoopEFTUncertainty                 /. output, 1.85778*^-05, eps];

Print["Number of passed tests: ", passedTests];
Print["Number of failed tests: ", failedTests];

Quit[failedTests];
