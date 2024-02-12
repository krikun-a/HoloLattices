#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Code:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*<<"IntermediateData/EquationCoefficients_back_lattice_1D.m";*)


gEOMsB


SetDirectory[NotebookDirectory[]];


<<"../../source/ImportMathematica/Pythonification.wl"


I gDyDzFT[[63]]^2+\[Mu]^2//Pythonify


(* ::Input::Initialization:: *)
Off[Part::pkspec1]
Off[Part::partd]


(* ::Input::Initialization:: *)
orderDiv={{0,0},{1,0},{0,1},{2,0},{0,2},{1,1}};
ne=6;


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]<>"Equations/"]*)


(* ::Subsubsection:: *)
(*Convert to Python*)


(* ::Input:: *)
(*Do[Module[{pyExpr},*)
(*pyExpr=(Pythonify/@Flatten[names[[1]]])/.{"gzi"->"z","gxi"->"z","gfb"->"f","gfi"->"f","gft"->"f"};*)
(*Export[names[[2]],StringRiffle[pyExpr,{"","\n"~~"",""}],"Text"]*)
(*],*)
(*{names,({{CoefEOMsB,"Coefs_1_0.txt"},{CoefEOMsT,"Coefs_1_1.txt"},{CoefEOMsI,"Coefs_I.txt"}})/.{gzi->z,gxi->z}}]*)
(*Do[Module[{pyExpr},*)
(*pyExpr=(Pythonify/@Flatten[names[[1]]])/.{gzi->z,gxi->z,gfb->f,gfi->f,gft->f};*)
(*Export[names[[2]],StringRiffle[pyExpr,{"","\n"~~"",""}],"Text"]*)
(*],*)
(*{names,({{gEOMsB,"EOMs_1_0.txt"},{gEOMsT,"EOMs_1_1.txt"},{gEOMsI,"EOMs_I.txt"}})}]*)


(* ::Input:: *)
(*Directory[]*)


gEOMsI//Variables
