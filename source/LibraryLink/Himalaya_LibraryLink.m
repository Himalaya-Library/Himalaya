BeginPackage["Himalaya`"];

InitializeHimalaya::usage = "Loads the Himalaya LibraryLink.  After
 InitializeHimalaya[] has been called, the function
 HimalayaCalculateDMh3L[] is available.

Usage:

  InitializeHimalaya[libraryLink];

Arguments:

 - libraryLink - The path to the Himalaya LibraryLink, 'Himalaya_LibraryLink.so'

";

HimalayaCalculateDMh3L::usage = "Runs Himalaya, given a set of input
parameters, and returns a list with all calculated loop corrections.

Usage:

  output = HimalayaCalculateDMh3L[settings -> {...}, parameters -> {...}];

Example:

MS = 2000;
TB = 20;
Xt = Sqrt[6]*MS;

output = HimalayaCalculateDMh3L[
    settings -> {
        bottom -> False,
        verbose -> True
    },
    parameters -> {
        scale -> MS,
        mu -> MS,
        g3 -> 0.1184,
        vd -> 246*Cos[ArcTan[TB]],
        vu -> 246*Sin[ArcTan[TB]],
        mq2 -> MS^2 IdentityMatrix[3],
        md2 -> MS^2 IdentityMatrix[3],
        mu2 -> MS^2 IdentityMatrix[3],
        Ab -> 0,
        At -> Xt + MS/TB,
        MA -> MS,
        MG -> MS,
        MW -> 80.384,
        MZ -> 91.1876,
        Mt -> 160,
        Mb -> 2.4
    }
];
";

{ parameters, settings };

Begin["`Private`"];

InitializeHimalaya[libName_String] := (
       HimalayaCalculateDMh3LLibInterface =
          LibraryFunctionLoad[libName, "HimalayaCalculateDMh3L", LinkObject, LinkObject];
    );

himalayaDefaultSettings = {
    bottom -> False,
    verbose -> True
};

himalayaDefaultParameters = {
        scale -> 0,
        mu -> 0,
        g3 -> 0,
        vd -> 0,
        vu -> 0,
        mq2 -> {{0,0,0},{0,0,0},{0,0,0}},
        md2 -> {{0,0,0},{0,0,0},{0,0,0}},
        mu2 -> {{0,0,0},{0,0,0},{0,0,0}},
        Ab -> 0,
        At -> 0,
        MA -> 0,
        MG -> 0,
        MW -> 0,
        MZ -> 0,
        Mt -> 0,
        Mb -> 0,
        MSt -> {-1, -1}, (* invalid input *)
        MSb -> {-1, -1}, (* invalid input *)
        s2t -> -2,       (* invalid input *)
        s2b -> -2        (* invalid input *)
};

Options[HimalayaCalculateDMh3L] = {
    Sequence @@ himalayaDefaultSettings,
    Sequence @@ himalayaDefaultParameters
};

Himalaya::nonum = "Error: `1` is not a numeric input value!";
Himalaya::error = "`1`";
Himalaya::info  = "`1`";
HimalayaInfoMessage[s_]  := Message[Himalaya::info, s];
HimalayaErrorMessage[s_] := Message[Himalaya::error, s];

HimalayaNumericQ[a_?NumericQ] := N[a];
HimalayaNumericQ[a_] := (Message[Himalaya::nonum, a]; Abort[]);

HimalayaCalculateDMh3L[a___, (settings | parameters) -> s_List, r___] :=
    HimalayaCalculateDMh3L[a, Sequence @@ s, r];

HimalayaCalculateDMh3L[OptionsPattern[]] :=
    HimalayaCalculateDMh3LLibInterface[
        HimalayaNumericQ /@ {
            (* settings *)
            Boole[OptionValue[bottom]],
            Boole[OptionValue[verbose]],
            (* parameters *)
            OptionValue[scale],
            OptionValue[mu],
            OptionValue[g3],
            OptionValue[vd],
            OptionValue[vu],
            OptionValue[mq2][[1,1]],
            OptionValue[mq2][[1,2]],
            OptionValue[mq2][[1,3]],
            OptionValue[mq2][[2,1]],
            OptionValue[mq2][[2,2]],
            OptionValue[mq2][[2,3]],
            OptionValue[mq2][[3,1]],
            OptionValue[mq2][[3,2]],
            OptionValue[mq2][[3,3]],
            OptionValue[md2][[1,1]],
            OptionValue[md2][[1,2]],
            OptionValue[md2][[1,3]],
            OptionValue[md2][[2,1]],
            OptionValue[md2][[2,2]],
            OptionValue[md2][[2,3]],
            OptionValue[md2][[3,1]],
            OptionValue[md2][[3,2]],
            OptionValue[md2][[3,3]],
            OptionValue[mu2][[1,1]],
            OptionValue[mu2][[1,2]],
            OptionValue[mu2][[1,3]],
            OptionValue[mu2][[2,1]],
            OptionValue[mu2][[2,2]],
            OptionValue[mu2][[2,3]],
            OptionValue[mu2][[3,1]],
            OptionValue[mu2][[3,2]],
            OptionValue[mu2][[3,3]],
            OptionValue[Ab],
            OptionValue[At],
            OptionValue[MA],
            OptionValue[MG],
            OptionValue[MW],
            OptionValue[MZ],
            OptionValue[Mt],
            OptionValue[Mb],
            OptionValue[MSt][[1]],
            OptionValue[MSt][[2]],
            OptionValue[MSb][[1]],
            OptionValue[MSb][[2]],
            OptionValue[s2t],
            OptionValue[s2b]
        }
    ];

End[];

EndPackage[];
