Get[FileNameJoin[{"..", "source", "LibraryLink", "Himalaya_LibraryLink.m"}]];

InitializeHimalaya["Himalaya_LibraryLink.so"];

MS = 2000;
TB = 20;
Xt = Sqrt[6] MS;

result = HimalayaCalculateDMh3L[
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
]

Print[result]
