(* ::Package:: *)

NotebookDirectory[]



(* ::Input::Initialization:: *)
Off[Part::pkspec1]
Off[Part::partd]


(* ::Input::Initialization:: *)
orderDiv={{0,0},{1,0},{0,1},{2,0},{0,2},{1,1}};
ne=6;


(* ::Input::Initialization:: *)
SubDerivArray={gFT[[i_]]:>gft[[i,1]],gDxFT[[i_]]:>gft[[i,2]],gDzFT[[i_]]:>gft[[i,3]],gD2xFT[[i_]]:>gft[[i,4]],gD2zFT[[i_]]:>gft[[i,5]],gDxDzFT[[i_]]:>gft[[i,6]],
gFI[[i_]]:>gfi[[i,1]],gDxFI[[i_]]:>gfi[[i,2]],gDzFI[[i_]]:>gfi[[i,3]],gD2xFI[[i_]]:>gfi[[i,4]],gD2zFI[[i_]]:>gfi[[i,5]],gDxDzFI[[i_]]:>gfi[[i,6]],
gFB[[i_]]:>gfb[[i,1]],gDxFB[[i_]]:>gfb[[i,2]],gDzFB[[i_]]:>gfb[[i,3]],gD2xFB[[i_]]:>gfb[[i,4]],gD2zFB[[i_]]:>gfb[[i,5]],gDxDzFB[[i_]]:>gfb[[i,6]],zeroI->0,zeroB->0,zeroT->0,gxI->gxi,gzI->gzi,gyI->gyi, dgDzFB[[i_]]:>dgfb[[i,3]],dgFB[[i_]]:>dgfb[[i,1]],dgD2zFB[[i_]]:>dgfb[[i,5]],dgD3zFB[[i_]]:>dgfb[[i,7]]}


(* ::Input::Initialization:: *)
Pythonify[expr0_]:=Module[{AdjustGreek, AdjustPart1, AdjustPart2,a1,a2,a3,a4,expr=expr0},
AdjustGreek={\[Mu]1->mu1,\[Mu]->mu, \[Omega]->omega,Sech[i_]:>1/cosh[i], Coth[i_]:>1/tanh[i],a_Complex:>Re[a]+HoldForm[Evaluate[Im[a]]]J};
AdjustPart1={"Part(":>"",","~~i:DigitCharacter..~~","~~j:DigitCharacter..~~")":>"["~~ToString[ToExpression[i]-1]~~","~~ToString[ToExpression[j]-1]~~"]"};
AdjustPart2={","~~i:DigitCharacter..~~")":>"["~~ToString[ToExpression[i]-1]~~"]"};
a1=expr/.AdjustGreek/.SubDerivArray;
a2=ToString[FortranForm[a1]];
a3=StringReplace[a2,AdjustPart1];
a3=StringReplace[a3,AdjustPart2];
a4=StringReplace[ToLowerCase[a3],{"j":>"1.J"}];
a4
]


(* ::Input::Initialization:: *)
Pythonify2D[expr0_]:=Module[{AdjustGreek, AdjustPart1, AdjustPart2,a1,a2,a3,a4,expr=expr0},
AdjustGreek={\[Mu]1->mu1,\[Mu]->mu,\[Omega]->omega, Sech[i_]:>1/cosh[i], Coth[i_]:>1/tanh[i],a_Complex:>Re[a]+HoldForm[Evaluate[Im[a]]]J};
AdjustPart1={"Part(":>"",","~~i:DigitCharacter..~~","~~j:DigitCharacter..~~")":>"["~~ToString[ToExpression[i]-1]~~","~~ToString[ToExpression[j]-1]~~"]"};
AdjustPart2={","~~i:DigitCharacter..~~")":>"["~~ToString[ToExpression[i]-1]~~"]"};
a1=expr/.AdjustGreek/.SubDerivArray2D;
a2=ToString[FortranForm[a1]];
a3=StringReplace[a2,AdjustPart1];
a3=StringReplace[a3,AdjustPart2];
a4=StringReplace[ToLowerCase[a3],{"j":>"1.J"}];
a4
]


Print["Pythonification loaded"]
