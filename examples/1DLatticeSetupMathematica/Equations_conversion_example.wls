#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Code:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*<<"IntermediateData/EquationCoefficients_back_lattice_1D.m";*)


gEOMsB


SetDirectory[NotebookDirectory[]];


NotebookDirectory[]


<<"../../source/ImportMathematica/Pythonification.wl"


I gDyDzFT[[63]]^2+\[Mu]^2//Pythonify


(* ::Input::Initialization:: *)
Off[Part::pkspec1]
Off[Part::partd]


orderDiv


ne


(* ::Input::Initialization:: *)
orderDiv={{0,0},{1,0},{0,1},{2,0},{0,2},{1,1}};
ne=6;


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]]*)
(*<<"/home/krikun/Bota/Works/2DLattice/EquationCoefficients_back_lattice_1D.m";*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]<>"IntermediateData/"]*)


(* ::Subsubsection:: *)
(*Convert to Python*)


(* ::Input:: *)
(*Do[Module[{pyExpr},*)
(*pyExpr=Pythonify/@Flatten[names[[1]]];*)
(*Export[names[[2]],StringRiffle[pyExpr,{"","\n"~~"",""}],"Text"]*)
(*],*)
(*{names,{{CoefEOMsB,"CoefsB_py.txt"},{CoefEOMsT,"CoefsT_py.txt"},{CoefEOMsI,"CoefsI_py.txt"}}}]*)
(*Do[Module[{pyExpr},*)
(*pyExpr=Pythonify/@Flatten[names[[1]]];*)
(*Export[names[[2]],StringRiffle[pyExpr,{"","\n"~~"",""}],"Text"]*)
(*],*)
(*{names,{{gEOMsB,"EOMsB_py.txt"},{gEOMsT,"EOMsT_py.txt"},{gEOMsI,"EOMsI_py.txt"}}}]*)


(* ::Input:: *)
(*Directory[]*)
