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
(*<<"EDCRGTCcode.m";*)
(*Off[General::partw];*)
(*Off[General::partd];*)
(*Off[General::pspec];*)


(* ::Input:: *)
(**)


(* ::Section:: *)
(*Compute the equations of motion*)


(* ::Input:: *)
(*Solve[T==(12-\[Mu]^2)/(16 \[Pi] \[Mu]),\[Mu]]*)


(* ::Subsection:: *)
(*Define coordinates, tetrad and interval*)


(* ::Input:: *)
(*Coord={t,x,y,z};*)
(*dCoord=d/@Coord;*)
(*Ndim=Length[Coord];*)


(* ::Input:: *)
(*$Assumptions=g0[z]>0&&f[z]>0&&P[z]>0&&Qtt[x,z]>0&&Qzz[x,z]>0&&Q[x,z]>0&&Qxz[x,z]>0&&M[x,z]>0&&Qty[x,z]>0&&R0[x,z]>0*)


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


(* ::Input:: *)
(*gRef=g/.{Qtt->(1&),Qzz->(1&),Q->(1&),M->(1&),Qxz->(0&),Qyz->(0&),Qty->(0&),Qtx->(0&),Qtz->(0&),R0->(0&)}*)
(**)


(* ::Subsection:: *)
(*Define a frame*)


(* ::Text:: *)
(*coframe has indices coframe[[a,mu]]=e^a_\mu*)
(*frame has indices frame[[a,mu]]=e_a^mu*)


(* ::Input::Initialization:: *)
\[Eta]={{-((f[z]Qtt[x,z])/g0[z]^2),0,0,0},{0,Q[x,z]/g0[z]^2 M[x,z],0,0},{0,0,Q[x,z]/g0[z]^2 M[x,z]^-1,0},{0,0,0,(g0'[z]^2 Qzz[x,z])/(f[z]g0[z]^2)}};


(* ::Input:: *)
(*coframe={{1, 0 Lx Qtx[x,z], 0 Ly Qty[x,z],0 Qtz[x,z]},{0,Lx,0,Qxz[x,z]},{0,0,Ly, 0Qyz[x,z]},{0,0,0,1}};*)
(*frame=Transpose[Inverse[coframe]];*)
(*gMUU = Transpose[frame] . \[Eta] . frame;*)
(*gMdd=Transpose[coframe] . \[Eta] . coframe;*)
(*MatrixForm[%]//Simplify//PowerExpand*)


(* ::Input:: *)
(*g-gMdd//Simplify*)


(* ::Subsection:: *)
(*Definer reference metric in the frame*)


(* ::Input:: *)
(*\[Eta]Ref=frame . gRef . Transpose[frame];*)
(*gMddRef=Transpose[coframe] . \[Eta]Ref . coframe;*)
(*MatrixForm[%]//Simplify//PowerExpand*)


(* ::Input:: *)
(*\[Eta]Ref//MatrixForm//Simplify*)


(* ::Subsection:: *)
(*Compute Christoffel for reference metric*)


(* ::Input:: *)
(*Variables[\[Eta]]*)


(* ::Input:: *)
(*Clear$dx*)
(*eIN=coframe . dCoord//Simplify*)


(* ::Input:: *)
(*Variables[eIN]*)


(* ::Input:: *)
(*Assuming[d[Lx]==0&&d[Ly]==0,RGtensors[\[Eta]Ref,Coord,eIN,{1,0,0}]];*)
(*GUddRef=GUdd;*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]<>"Data/"]*)


(* ::Input:: *)
(*GUddRef>>"GUddRef_Lattice1D_tetrad_Donos.dat";*)
(**)


(* ::Subsection:: *)
(*RG computation*)


(* ::Input:: *)
(*Clear$dx*)
(*eIN=coframe . dCoord//Simplify*)


(* ::Input:: *)
(*Assuming[d[Lx]==0&&d[Ly]==0,RGtensors[\[Eta],Coord,eIN,{1,0,0}]];*)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*{gUU,GUdd,gdd,Rdd}>>"Data/Grav_data_1D_lattice_tetrad_Donos.dat";*)


(* ::Subsection:: *)
(*Define DeTurck vector*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*\[Xi]U=multiDot [gUU,(GUdd-GUddRef),{1,2},{2,3}];*)
(*\[Xi]d =Lower[\[Xi]U,1];*)


(* ::Input:: *)
(*\[Xi]Norm=\[Xi]d . \[Xi]U//Refine;*)
(**)


(* ::Input:: *)
(*divXi = (1/2 (Transpose[#]+#)&)@covD[\[Xi]d];*)
(**)


(* ::Subsection:: *)
(*Definitions of fields and potential*)


(* ::Input:: *)
(*?Sech*)


(* ::Input:: *)
(*\[Tau]=(Sech[Sqrt[3]#]&);*)
(*V=(-6 Cosh[Sqrt[1/3]#]&);*)
(*W=(V[#]+6&);*)
(*\[Kappa]=1/2;*)
(*\[Theta]=(c1/(6Sqrt[2]) Tanh[Sqrt[3]#]&);*)
(*Clear[\[Tau],V,W,\[Kappa],\[Theta]]*)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*$Assumptions=$Assumptions&&\[Psi][x,z]>0&&at[x,z]>0&&ax[x,z]>0&&ay[x,z]>0&&az[x,z]>0*)


(* ::Input:: *)
(*Ad =frame . {1/4 \[Mu] g0'[z]^2 at[x,z], 0,0,0};*)


(* ::Input:: *)
(*Fdd=(Transpose[#]-#&)@ covD[Ad];*)
(**)


(* ::Input:: *)
(*Fdd//Simplify//MatrixForm*)
(**)


(* ::Input:: *)
(*PsiA=covDiv[gUU . Ad,{1}] + \[Xi]U . Ad;*)
(*dPsiA=covD[PsiA];*)


(* ::Input:: *)
(*Coefficient[EOMAH[[4]],Derivative[0,0,2][az][x,y,z]]*)


(* ::Subsection:: *)
(*Equations of motion*)


(* ::Input:: *)
(*Gdd=Rdd+3gdd  - 1/2 (multiDot[Fdd,Raise[Fdd,2],{2,2}] - gdd/4 multiDot[Fdd,Raise[Fdd,{1,2}],{1,1},{2,2}])//Refine;*)


(* ::Input:: *)
(*EOMA=covDiv[Raise[Fdd,1],{1}]//Refine;*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*{Gdd,EOMA}>>"Data/Lattice_1D_raw_EOMs_tetrad_Donos.dat";*)


(* ::Subsection:: *)
(*DeTurck equation*)


(* ::Input:: *)
(*GHdd = Gdd - divXi//Refine;*)


(* ::Input:: *)
(*EOMAH=EOMA+1dPsiA//Refine;*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*{GHdd,EOMAH,divXi,dPsiA}>>"Data/Lattice_1D_gauge_EOMs_tetrad_Donos.dat";*)


(* ::Input:: *)
(*dPsiA[[1]]*)


(* ::Subsection:: *)
(*Check the nonzero entries in the equations*)


(* ::Input:: *)
(*Map[If[!NumberQ[#],1,0]&,GHdd,{2}]//MatrixForm*)


(* ::Input:: *)
(*Map[If[!NumberQ[#],1,0]&,EOMAH,{1}]//MatrixForm*)


(* ::Input:: *)
(*Map[If[!NumberQ[#],1,0]&,\[Xi]d,{1}]//MatrixForm*)
(**)


(* ::Input:: *)
(*trG=Contract[gUU . Gdd,{1,2}]//Refine;*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*{\[Xi]Norm,trG}>>"Data/XiNorm_TrG_1D_lattice_tetrad_Donos.dat";*)


(* ::Subsection:: *)
(*Normalize*)


(* ::Input:: *)
(*Variables[GHdd]*)


(* ::Input::Initialization:: *)
Bmodes={at[x,z],Qtt[x,z],M[x,z],Q[x,z],Qzz[x,z],Qxz[x,z]};


(* ::Input:: *)
(*Length[Bmodes]*)


(* ::Input::Initialization:: *)
DBmodes=Flatten[Through[{D[#,x]&,D[#,z]&,D[#,x,x]&,D[#,z,z]&,D[#,x,z]&}[Bmodes]]];


(* ::Input:: *)
(*DBmodes[[6*3+1;;6*4]]*)


(* ::Input:: *)
(*BeqsRaw={EOMAH[[1]],GHdd[[1,1]],GHdd[[2,2]]-M[x,z]^2 GHdd[[3,3]],1/2 (GHdd[[2,2]]+M[x,z]^2 GHdd[[3,3]]),GHdd[[4,4]],1/M[x,z] GHdd[[4,2]]};*)


(* ::Input:: *)
(*Transpose[Coefficient[BeqsRaw,#]&/@DBmodes[[6*3+1;;6*4]]]//Simplify//MatrixForm*)
(**)


(* ::Input:: *)
(*BeqsSimplified=ConstantArray[0,6]*)


(* ::Input:: *)
(*Table[time=Timing[BeqsSimplified[[j]]=Collect[BeqsRaw[[j]]/Coefficient[BeqsRaw[[j]],DBmodes[[6*3+j]]],Join[DBmodes,{Bmodes[[j]]}],(Simplify[#,TimeConstraint->{30,30}]&)];*)
(*];*)
(*Print["j="<>ToString[j]<>" t="<>ToString[time[[1]]]];*)
(*,{j,1,6}]*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*BeqsSimplified>>"Data/Lattice1D_EOMs_simplified_tetrad_Donos.dat";*)


(* ::Input:: *)
(*DBmodes[[8*4+1;;8*5]]*)


(* ::Input:: *)
(*Transpose[Coefficient[BeqsSimplified,#]&/@DBmodes[[15*3+1;;15*4]]]//Simplify//MatrixForm*)


(* ::Input::Initialization:: *)
SubF={at->( F[1][#1,#2]&),Qtt->F[2],M->F[3],Q->F[4],Qzz->F[5],Qxz->F[6]};


(* ::Subsection:: *)
(*Rescale coordinates (already done in the ansatz)*)


(* ::Text:: *)
(*c.f. Donos 1512.06861 (5.2)*)


(* ::Input::Initialization:: *)
SubFxz=F[i_]:>Fxz[[i]];

SubPotentials=MapThread[(#1->#2&),{{f,g0},{((1-#)^2 (1+(1-(1-#)^2)+(1-(1-#)^2)^2- (\[Mu]1^2 (1-(1-#)^2)^3)/4)&),(1-(1-#)^2&)}}];


(* ::Input:: *)
(*1/4 (\[Mu]1^2 (r-2)^3 r^3+4 (1+2r+3 r^2-4 r^3+r^4))-(1+(1-(1-r)^2)+(1-(1-r)^2)^2- (\[Mu]1^2 (1-(1-r)^2)^3)/4)//Simplify*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*BeqsSimplified[[6]]/.SubF/.SubFxz/.SubPotentials//PowerExpand;*)
(*Variables[%]*)


(* ::Input:: *)
(*EOMsI=BeqsSimplified/.SubF/.SubFxz;*)


(* ::Input::Initialization:: *)
SubLinx=Fxz[[i_]]:>(Fxz[[i]][#1,#2]+lin dFxz[[i]][#1,#2]&);



(* ::Input:: *)
(*dEOMsI=(D[EOMsI/.SubLinx,lin])/.{lin->0};*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*EOMsI>>"Data/EOMsI_lattice_1D_tetrad_Donos.dat";*)


(* ::Input:: *)
(*DumpSave["Data/dEOMsI_lattice_1D_tetrad_Donos.m",dEOMsI];*)


(* ::Input:: *)
(*Coefficient[EOMsI[[1]],Derivative[0,2][Fxz[[1]]][x,z]]*)


(* ::Subsection:: *)
(*Boundary BC*)


(* ::Input:: *)
(*Bmodes*)


(* ::Input:: *)
(*EOMsB=Table[Fxz[[i]][x,z]-If[i>=2&&i<=5,1,0],{i,1,6}];*)
(*EOMsB[[1]]-=(1 +Ax Cos[NperiodsX 2 Pi x+Pi phaseX]);*)
(*EOMsB;*)


(* ::Input:: *)
(*dEOMsB=(D[EOMsB/.SubLinx,lin])/.{lin->0};*)


(* ::Input:: *)
(*EOMsB*)


(* ::Subsection:: *)
(*Horizon BC*)


(* ::Input:: *)
(*Bmodes*)


(* ::Input:: *)
(*Bmodes[[5]]*)


(* ::Input:: *)
(*EOMsT=Table[Derivative[0,1][Fxz[[i]]][x,z],{i,1,6}];*)
(*EOMsT[[5]]=Fxz[[5]][x,z]-Fxz[[2]][x,z];*)
(*EOMsT*)


(* ::Input:: *)
(*dEOMsT=(D[EOMsT/.SubLinx,lin])/.{lin->0};*)
(**)


(* ::Subsection:: *)
(*Coefficient arrays*)


(* ::Input:: *)
(*orderDiv={{0,0},{1,0},{0,1},{2,0},{0,2},{1,1}};*)
(*Length[orderDiv]*)
(*ne=6;*)


(* ::Input:: *)
(*SubFI={z -> gzI, x -> gxI,*)
(*Fxz[[i_]][__] :> gFI[[i]],*)
(*Derivative[1, 0][Fxz[[i_]]][__] :> gDxFI[[i]],*)
(*Derivative[0, 1][Fxz[[i_]]][__] :> gDzFI[[i]], *)
(*Derivative[2, 0][Fxz[[i_]]][__] :> gD2xFI[[i]],*)
(*Derivative[0, 2][Fxz[[i_]]][__] :> gD2zFI[[i]],*)
(*Derivative[1, 1][Fxz[[i_]]][__] :> gDxDzFI[[i]]};*)
(*SubFB={z -> gzB, x -> gxB,*)
(*Fxz[[i_]][__] :> gFB[[i]],*)
(*Derivative[1, 0][Fxz[[i_]]][__] :> gDxFB[[i]],*)
(*Derivative[0, 1][Fxz[[i_]]][__] :> gDzFB[[i]], *)
(*Derivative[2, 0][Fxz[[i_]]][__] :> gD2xFB[[i]],*)
(*Derivative[0, 2][Fxz[[i_]]][__] :> gD2zFB[[i]],*)
(*Derivative[1, 1][Fxz[[i_]]][__] :> gDxDzFB[[i]]};*)
(*SubFT={z -> gzT, x -> gxT,*)
(*Fxz[[i_]][__] :> gFT[[i]],*)
(*Derivative[1, 0][Fxz[[i_]]][__] :> gDxFT[[i]],*)
(*Derivative[0, 1][Fxz[[i_]]][__] :> gDzFT[[i]], *)
(*Derivative[2, 0][Fxz[[i_]]][__] :> gD2xFT[[i]],*)
(*Derivative[0, 2][Fxz[[i_]]][__] :> gD2zFT[[i]],*)
(*Derivative[1, 1][Fxz[[i_]]][__] :> gDxDzFT[[i]]};*)


(* ::Input:: *)
(*dEOMsI;*)
(*Map[Coefficient[#,Flatten[Array[Derivative[orderDiv[[#1,1]],orderDiv[[#1,2]]][dFxz[[#2]]][x,z]&,{6,ne}]]]&,%];*)
(*CoefEOMsI=zeroI+Simplify[%,TimeConstraint->10] /. SubFI/.SubPotentials;*)


(* ::Input:: *)
(*EOMsI/. SubFI/.SubPotentials;*)
(*gEOMsI=%;*)


(* ::Input:: *)
(*dEOMsT;*)
(*Map[Coefficient[#,Flatten[Array[Derivative[orderDiv[[#1,1]],orderDiv[[#1,2]]][dFxz[[#2]]][x,z]&,{6,ne}]]]&,%];*)
(*CoefEOMsT=zeroT+Simplify[%,TimeConstraint->10] /. SubFT/.SubPotentials;*)


(* ::Input:: *)
(*EOMsT/. SubFT/.SubPotentials;*)
(*gEOMsT=%;*)


(* ::Input:: *)
(*dEOMsB;*)
(*Map[Coefficient[#,Flatten[Array[Derivative[orderDiv[[#1,1]],orderDiv[[#1,2]]][dFxz[[#2]]][x,z]&,{6,ne}]]]&,%];*)
(*CoefEOMsB=zeroB+Simplify[%,TimeConstraint->10] /. SubFB/.SubPotentials;*)
(*EOMsB/.SubFB/.SubPotentials;*)
(*gEOMsB=%;*)


(* ::Input:: *)
(*Variables[gEOMsI]*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*CoefEOMsI>>"Data/CoefEOMsI_1D_no_potentials_tetrad_Donos.dat";*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*<<"!rm EquationCoefficients_back_lattice_1D.m"*)
(*Save["EquationCoefficients_back_lattice_1D.m",{ne,orderDiv,gEOMsB,gEOMsI,gEOMsT,CoefEOMsB,CoefEOMsI,CoefEOMsT}]*)
(**)


(* ::Input:: *)
(*Variables[{gEOMsI,gEOMsB}]*)


(* ::Subsection:: *)
(*Control structures*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*{\[Xi]Norm,trG}=<<"Data/XiNorm_TrG_1D_lattice_tetrad_Donos.dat";*)


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
(*<<"!rm ControlExpressions_back_lattice_1D_dump.m";*)
(*DumpSave["ControlExpressions_back_lattice_1D_dump.m",{gNormXiI,gNormXiT,gNormXiB,gtrGI,gtrGT,gtrGB}];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*<<"!rm ControlExpressions_back_lattice_1D.m";*)
(*Save["ControlExpressions_back_lattice_1D.m",{gNormXiI,gNormXiT,gNormXiB,gtrGI,gtrGT,gtrGB}];*)
