#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Preamble*)


(* ::Subsection:: *)
(*Description*)


(* ::Text:: *)
(*This notebook demonstrates how one obtains the equations of motion and the appropriate boundary conditions for the holographic treatment of a 2D metal with ionic lattice in 1 direction. *)
(*The model consists of the Abelian vector field in presence of dynamical gravity in the planar AdS black hole background and with a nontrivial boundary condition, defining the explicit periodic lattice.*)
(*This setup is discussed in arXiv:1209.1098.*)
(*The particular realization with tetrads, which we use here, is inherited from 1512.06861.*)
(*Requirements: EDCRGTC*)


(* ::Subsection:: *)
(*Load EDCRGTC*)


(* ::Text:: *)
(*We will use EDCRGTC for evaluating the gravity equations of motion*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*<<"EDCRGTCcode.m"*)
(*Off[General::partw];*)
(*Off[General::partd];*)
(*Off[General::pspec];*)
(*On[Assert];*)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Compute the equations of motion*)


(* ::Input:: *)
(*Solve[T==(12-\[Mu]^2)/(16 \[Pi] \[Mu]),\[Mu]]*)


(* ::Subsection:: *)
(*Geometry setup*)


(* ::Subsubsection:: *)
(*Define coordinates, interval and metric*)


(* ::Input:: *)
(*Coord={t,x,y,z};*)
(*dCoord=d/@Coord;*)
(*Ndim=Length[Coord];*)


(* ::Input:: *)
(*$Assumptions=g0[z]>0&&f[z]>0&&P[z]>0&&Qtt[x,z]>0&&Qzz[x,z]>0&&Q[x,z]>0&&Qxz[x,z]>0&&M[x,z]>0&&Qty[x,z]>0&&R0[x,z]>0*)


(* ::Text:: *)
(*The metric ansatz is a reduced version of 1512.06861.*)


(* ::Input:: *)
(*\[Eta]T=d[t];*)
(*\[Eta]X=Lx d[x] + Qxz[x,z] d[z] ;*)
(*\[Eta]Y=Ly d[y];*)


(* ::Input:: *)
(*interval=1/g0[z]^2 (-f[z]Qtt[x,z] (\[Eta]T)^2 + (g0'[z]^2 Qzz[x,z] d[z]^2)/f[z]+Q[x,z](M[x,z](\[Eta]X)^2+M[x,z]^(-1) (\[Eta]Y)^2) );*)
(*ds=Collect[interval,dCoord,Simplify];*)
(*Map[Coefficient[ds,#]&,  TensorProduct[dCoord,dCoord],{2}];*)
(*g=Table[If[i==j,%[[i,j]],%[[i,j]]/2],{i,1,Ndim},{j,1,Ndim}];*)


(* ::Input:: *)
(*g//MatrixForm*)


(* ::Text:: *)
(*The reference metric for the DeTurck trick*)


(* ::Input:: *)
(*gRef=g/.{Qtt->(1&),Qzz->(1&),Q->(1&),M->(1&),Qxz->(0&)};*)


(* ::Subsubsection:: *)
(*Define a frame*)


(* ::Text:: *)
(*Coframe has the indices coframe[[a,mu]]=e^a_\mu*)
(*Frame has indices frame[[a,mu]]=e_a^mu*)


(* ::Input:: *)
(*\[Eta]={{-((f[z]Qtt[x,z])/g0[z]^2),0,0,0},{0,Q[x,z]/g0[z]^2 M[x,z],0,0},{0,0,Q[x,z]/g0[z]^2 M[x,z]^-1,0},{0,0,0,(g0'[z]^2 Qzz[x,z])/(f[z]g0[z]^2)}};*)


(* ::Input:: *)
(*coframe={{1,0, 0,0},{0,Lx,0,Qxz[x,z]},{0,0,Ly, 0},{0,0,0,1}};*)
(*frame=Transpose[Inverse[coframe]];*)
(*gMUU = Transpose[frame] . \[Eta] . frame;*)
(*gMdd=Transpose[coframe] . \[Eta] . coframe;*)


(* ::Text:: *)
(*Check that the product of tetrads equals the metric*)


(* ::Input:: *)
(*Assert[Simplify[Total[Flatten[g-gMdd]^2]]==0]*)


(* ::Subsubsection:: *)
(*Define the reference metric in the frame*)


(* ::Input:: *)
(*\[Eta]Ref=frame . gRef . Transpose[frame];*)
(*gMddRef=Transpose[coframe] . \[Eta]Ref . coframe;*)


(* ::Input:: *)
(*\[Eta]Ref//MatrixForm//Simplify;*)


(* ::Subsubsection:: *)
(*Compute Christoffel for reference metric*)


(* ::Input:: *)
(*Clear$dx;*)
(*eIN=coframe . dCoord//Simplify;*)


(* ::Input:: *)
(*Assuming[d[Lx]==0&&d[Ly]==0,RGtensors[\[Eta]Ref,Coord,eIN,{1,0,0}]];*)
(*GUddRef=GUdd;*)


(* ::Subsubsection:: *)
(*Compute Christoffel for the ansatz metric*)


(* ::Input:: *)
(*Clear$dx*)
(*eIN=coframe . dCoord//Simplify;*)


(* ::Input:: *)
(*Assuming[d[Lx]==0&&d[Ly]==0,RGtensors[\[Eta],Coord,eIN,{1,0,0}]];*)
(**)


(* ::Subsubsection:: *)
(*Define DeTurck vector*)


(* ::Input:: *)
(*\[Xi]U=multiDot [gUU,(GUdd-GUddRef),{1,2},{2,3}];*)
(*\[Xi]d =Lower[\[Xi]U,1];*)


(* ::Input:: *)
(*\[Xi]Norm=\[Xi]d . \[Xi]U//Refine;*)


(* ::Input:: *)
(*divXi = (1/2 (Transpose[#]+#)&)@covD[\[Xi]d];*)


(* ::Subsection:: *)
(*Equations of motion*)


(* ::Subsubsection:: *)
(*Definitions of fields and potentials*)


(* ::Input:: *)
(*$Assumptions=$Assumptions&&\[Psi][x,z]>0&&at[x,z]>0&&ax[x,z]>0&&ay[x,z]>0&&az[x,z]>0*)


(* ::Input:: *)
(*Ad =frame . {1/4 \[Mu] g0'[z]^2 at[x,z], 0,0,0};*)


(* ::Input:: *)
(*Fdd=(Transpose[#]-#&)@ covD[Ad];*)


(* ::Text:: *)
(*Gauge fixing term for the Maxwell field*)


(* ::Input:: *)
(*PsiA=covDiv[gUU . Ad,{1}] + \[Xi]U . Ad;*)
(*dPsiA=covD[PsiA];*)


(* ::Subsubsection:: *)
(*Einstein equation*)


(* ::Input:: *)
(*Gdd=Rdd+3gdd  - 1/2 (multiDot[Fdd,Raise[Fdd,2],{2,2}] - gdd/4 multiDot[Fdd,Raise[Fdd,{1,2}],{1,1},{2,2}])//Refine;*)


(* ::Text:: *)
(*DeTurck equation*)


(* ::Input:: *)
(*GHdd = Gdd - divXi//Refine;*)


(* ::Subsubsection:: *)
(*Maxwell equation*)


(* ::Input:: *)
(*EOMA=covDiv[Raise[Fdd,1],{1}]//Refine;*)


(* ::Text:: *)
(*Gauge fixed Maxwell equation*)


(* ::Input:: *)
(*EOMAH=EOMA+dPsiA//Refine;*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*{GHdd,EOMAH,divXi,dPsiA}>>"IntermediateData/Lattice_1D_gauge_EOMs.dat";*)


(* ::Subsubsection:: *)
(*Check the nonzero entries in the equations*)


(* ::Input:: *)
(*Map[If[!NumberQ[#],1,0]&,GHdd,{2}]//MatrixForm*)


(* ::Input:: *)
(*Map[If[!NumberQ[#],1,0]&,EOMAH,{1}]//MatrixForm*)


(* ::Input:: *)
(*Map[If[!NumberQ[#],1,0]&,\[Xi]d,{1}]//MatrixForm*)


(* ::Subsubsection:: *)
(*Define the control expressions*)


(* ::Input:: *)
(*trG=Contract[gUU . Gdd,{1,2}]//Refine;*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*{\[Xi]Norm,trG}>>"IntermediateData/XiNorm_TrG.dat";*)


(* ::Subsection:: *)
(*Normalize*)


(* ::Input:: *)
(*Variables[GHdd]*)


(* ::Text:: *)
(*List and count the ansatz functions we are solving for*)


(* ::Input::Initialization:: *)
Bmodes={at[x,z],Qtt[x,z],M[x,z],Q[x,z],Qzz[x,z],Qxz[x,z]};
ne=Length[Bmodes]


(* ::Input::Initialization:: *)
DBmodes=Flatten[Through[{D[#,x]&,D[#,z]&,D[#,x,x]&,D[#,z,z]&,D[#,x,z]&}[Bmodes]]];


(* ::Text:: *)
(*Diagonalize the equations*)


(* ::Input:: *)
(*BeqsRaw={EOMAH[[1]],GHdd[[1,1]],GHdd[[2,2]]-M[x,z]^2 GHdd[[3,3]],1/2 (GHdd[[2,2]]+M[x,z]^2 GHdd[[3,3]]),GHdd[[4,4]],1/M[x,z] GHdd[[4,2]]};*)
(*Transpose[Coefficient[BeqsRaw,#]&/@DBmodes[[6*3+1;;6*4]]]//Simplify//MatrixForm*)


(* ::Text:: *)
(*Normalize the leading derivative terms*)


(* ::Input:: *)
(*BeqsSimplified=ConstantArray[0,ne];*)
(*Table[BeqsSimplified[[j]]=Collect[BeqsRaw[[j]]/Coefficient[BeqsRaw[[j]],DBmodes[[ne*3+j]]],Join[DBmodes,{Bmodes[[j]]}],(Simplify[#,TimeConstraint->{30,30}]&)];*)
(*,{j,1,ne}];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*BeqsSimplified>>"IntermediateData/Lattice1D_EOMs_normalized.dat"*)


SetDirectory[NotebookDirectory[]];
BeqsSimplified=<<"IntermediateData/Lattice1D_EOMs_normalized.dat";


(* ::Subsection:: *)
(*Rename the functions to solve for*)


(* ::Input::Initialization:: *)
SubF={at->( F[1][#1,#2]&),Qtt->F[2],M->F[3],Q->F[4],Qzz->F[5],Qxz->F[6]};


(* ::Subsection:: *)
(*Rescale coordinates (already done in the ansatz)*)


(* ::Text:: *)
(*c.f. Donos 1512.06861 (5.2)*)


(* ::Input::Initialization:: *)
SubFxz=F[i_]:>Fxz[[i]];
SubPotentials=MapThread[(#1->#2&),{{f,g0},{((1-#)^2 (1+(1-(1-#)^2)+(1-(1-#)^2)^2- (\[Mu]1^2 (1-(1-#)^2)^3)/4)&),(1-(1-#)^2&)}}];


(* ::Subsection:: *)
(*Internal equations and derivatives*)


(* ::Input:: *)
(*EOMsI=BeqsSimplified/.SubF/.SubFxz;*)


(* ::Input:: *)
(*SubLinx=Fxz[[i_]]:>(Fxz[[i]][#1,#2]+lin dFxz[[i]][#1,#2]&);*)


(* ::Input:: *)
(*dEOMsI=(D[EOMsI/.SubLinx,lin])/.{lin->0};*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*EOMsI>>"IntermediateData/EOMsI.dat";*)


(* ::Input:: *)
(*DumpSave["IntermediateData/dEOMsI.m",dEOMsI];*)


SetDirectory[NotebookDirectory[]]
EOMsI=<<"IntermediateData/EOMsI.dat";
<<"IntermediateData/dEOMsI.m";


dEOMsI//Variables;


(* ::Subsection:: *)
(*Boundary BC*)


SubF


(* ::Text:: *)
(*Unity for all the metric fields except off-diagonal Qxz. This sets asymptotic AdS.*)


(* ::Input:: *)
(*EOMsB=Table[Fxz[[i]][x,z]-If[i>=2&&i<=5,1,0],{i,1,ne}];*)


(* ::Text:: *)
(*Here we define the shape of the boundary lattice*)


(* ::Input:: *)
(*EOMsB[[1]]-=(1 +Ax Cos[NperiodsX 2 Pi x+Pi phaseX]);*)


(* ::Input:: *)
(*dEOMsB=(D[EOMsB/.SubLinx,lin])/.{lin->0};*)


(* ::Subsection:: *)
(*Horizon BC*)


(* ::Text:: *)
(*For the rescaled radial coordinate all derivatives are zero on the horizon*)


(* ::Input:: *)
(*EOMsT=Table[Derivative[0,1][Fxz[[i]]][x,z],{i,1,6}];*)


(* ::Text:: *)
(*Thermal equilibrium requires Q_tt = Q_zz*)


(* ::Input:: *)
(*EOMsT[[5]]=Fxz[[5]][x,z]-Fxz[[2]][x,z];*)


(* ::Input:: *)
(*dEOMsT=(D[EOMsT/.SubLinx,lin])/.{lin->0};*)


(* ::Subsection:: *)
(*Coefficient arrays*)


(* ::Input:: *)
(*orderDiv={{0,0},{1,0},{0,1},{2,0},{0,2},{1,1}};*)
(*ndiv=Length[orderDiv]*)


(* ::Input:: *)
(*SubFI={z -> gzI, x -> gxI,*)
(*Fxz[[i_]][__] :> gFI[[i]],*)
(*Derivative[1, 0][Fxz[[i_]]][__] :> gDxFI[[i]],*)
(*Derivative[0, 1][Fxz[[i_]]][__] :> gDzFI[[i]], *)
(*Derivative[2, 0][Fxz[[i_]]][__] :> gD2xFI[[i]],*)
(*Derivative[0, 2][Fxz[[i_]]][__] :> gD2zFI[[i]],*)
(*Derivative[1, 1][Fxz[[i_]]][__] :> gDxDzFI[[i]]};*)
(**)


(* ::Input:: *)
(*dEOMsI;*)
(*Map[Coefficient[#,Flatten[Array[Derivative[orderDiv[[#1,1]],orderDiv[[#1,2]]][dFxz[[#2]]][x,z]&,{ndiv,ne}]]]&,%];*)
(*CoefEOMsI=zeroI + Simplify[%,TimeConstraint->10] /. SubFI/.SubPotentials;*)


(* ::Input:: *)
(*EOMsI/. SubFI/.SubPotentials;*)
(*gEOMsI=%;*)


(* ::Input:: *)
(*dEOMsT;*)
(*Map[Coefficient[#,Flatten[Array[Derivative[orderDiv[[#1,1]],orderDiv[[#1,2]]][dFxz[[#2]]][x,z]&,{ndiv,ne}]]]&,%];*)
(*CoefEOMsT=zeroT+Simplify[%,TimeConstraint->10] /. SubFI/.SubPotentials;*)


(* ::Input:: *)
(*EOMsT/. SubFI/.SubPotentials;*)
(*gEOMsT=%;*)


(* ::Input:: *)
(*dEOMsB;*)
(*Map[Coefficient[#,Flatten[Array[Derivative[orderDiv[[#1,1]],orderDiv[[#1,2]]][dFxz[[#2]]][x,z]&,{ndiv,ne}]]]&,%];*)
(*CoefEOMsB=zeroB+Simplify[%,TimeConstraint->10] /. SubFI/.SubPotentials;*)
(*EOMsB/.SubFI/.SubPotentials;*)
(*gEOMsB=%;*)


(* ::Input:: *)
(*Variables[gEOMsB];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*<<"!rm IntermediateData/EquationCoefficients_back_lattice_1D.m"*)
(*Save["IntermediateData/EquationCoefficients_back_lattice_1D.m",{ne,orderDiv,gEOMsB,gEOMsI,gEOMsT,CoefEOMsB,CoefEOMsI,CoefEOMsT}]*)


(* ::Subsection:: *)
(*Control structures*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*{\[Xi]Norm,trG}=<<"IntermediateData/XiNorm_TrG.dat";*)


(* ::Input:: *)
(*XiNormI=\[Xi]Norm/z^2/.SubF/.SubFxz/.SubPotentials//Refine;*)
(*XiNormT=0;*)
(*XiNormB=0;*)
(*trGI=trG/z/.SubF/.SubFxz/.SubPotentials//Refine;*)
(*trGT=0;*)
(*trGB=0;*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*gNormXiI=XiNormI/.SubFI;*)
(*gNormXiT=XiNormT/.SubFT;*)
(*gNormXiB=XiNormB/.SubFB;*)
(**)


(* ::Input:: *)
(*gNormXiB=zeroB;*)
(*gtrGB=zeroB;*)
(*gNormXiT=zeroT;*)
(*gtrGT=zeroT;*)


(* ::Input:: *)
(*gtrGI=trGI/.SubFI;*)
(*gtrGT=trGT/.SubFT;*)
(*gtrGB=trGB/.SubFB;*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Variables[gEOMsB]*)


(* ::Input:: *)
(*Variables[gEOMsI]*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*<<"!rm IntermediateData/ControlExpressions_back_lattice_1D_dump.m";*)
(*DumpSave["IntermediateData/ControlExpressions_back_lattice_1D_dump.m",{gNormXiI,gNormXiT,gNormXiB,gtrGI,gtrGT,gtrGB}];*)
