(* ::Package:: *)

BeginPackage["fourdiff`"];


fourdiff::usage = "Differentiate periodic function wrt t";
mik::usage = "Multiply Fourier modes by -ik ";
mikold::usage = "Multiply Fourier modes by -ik ";
double::usage = "Double fourier modes for anti-aliasing";
doubleold::usage = "Double fourier modes for anti-aliasing";
halve::usage = "Halve number of fourier modes for anti-aliasing.";
halveold::usage = "Halve number of fourier modes for anti-aliasing.";


Begin["`Private`"];


fourdiff[f_,Delta_,Nt_]:=Module[{faux,f0,bigf},
faux=1/Sqrt[Nt]*Fourier[f];
faux=mik[faux,Delta,Nt];
faux=Sqrt[Nt]*InverseFourier[faux];
faux
];


mikold[fmodes_,Delta_,Nt_]:=Module[{dfmodes},
(*positive frequency*)
dfmodes=ConstantArray[0.,Nt];
Do[dfmodes[[k+1]]=-I*k*(2\[Pi])/Delta *fmodes[[k+1]],{k,0,IntegerPart[Nt/2-1]}];
(*high frequency cosine*)
dfmodes[[IntegerPart[Nt/2+1]]]=0;
(*negative frequency*)
Do[dfmodes[[k+1]]=fmodes[[k+1]]*(-I*(k-Nt)*(2\[Pi])/Delta),{k,IntegerPart[Nt/2+1],Nt-1}];
dfmodes
];


mik=Compile[{{fmodes,_Complex,1},{Delta,_Real},{Nt,_Integer}},Module[{dfmodes},
(*positive frequency*)
dfmodes=fmodes;
Do[dfmodes[[k+1]]=-I*k*(2\[Pi])/Delta *fmodes[[k+1]],{k,0,IntegerPart[Nt/2-1]}];
(*high frequency cosine*)
dfmodes[[IntegerPart[Nt/2+1]]]=0;
(*negative frequency*)
Do[dfmodes[[k+1]]=fmodes[[k+1]]*(-I*(k-Nt)*(2\[Pi])/Delta),{k,IntegerPart[Nt/2+1],Nt-1}];
dfmodes
]
,CompilationTarget->"C"];


doubleold[fmodes_,Nt_]:=Module[{faux},
	faux=ConstantArray[0.,2*Nt];
	(*positive frequencies*)
	Do[faux[[k+1]]=fmodes[[k+1]],{k,0,Nt/2-1}];
	(*split high-frequency cosine*)
	faux[[Nt/2+1]]=1/2*fmodes[[Nt/2+1]];
	faux[[3*Nt/2+1]]=1/2*fmodes[[Nt/2+1]];
	(*negative frequencies*)
	Do[faux[[k+Nt+1]]=fmodes[[k+1]],{k,Nt/2+1,Nt-1}];
	faux
];


double=Compile[{{fmodes,_Complex,1},{Nt,_Integer}}, Module[{faux},
	(*Initialize array, it is unimportant how it is filled*)
	faux=ConstantArray[0.+I*0.,2*Nt];
	(*positive frequencies*)
	Do[faux[[k+1]]=fmodes[[k+1]],{k,0,IntegerPart[Nt/2-1]}];
	(*split high-frequency cosine*)
	faux[[IntegerPart[Nt/2+1]]]=1./2.*fmodes[[IntegerPart[Nt/2+1]]];
	faux[[IntegerPart[3*Nt/2+1]]]=1./2.*fmodes[[IntegerPart[Nt/2+1]]];
	(*negative frequencies*)
	Do[faux[[k+Nt+1]]=fmodes[[k+1]],{k,IntegerPart[Nt/2+1],Nt-1}];
	faux
]
,CompilationTarget->"C"];



halveold[fmodes_,Nt_]:=Module[{faux},
	faux=ConstantArray[0.,Nt];
	(*positive frequencies*)
	Do[faux[[k+1]]=fmodes[[k+1]],{k,0,Nt/2-1}];
	(*sum high-frequency cosines*)
	faux[[Nt/2+1]]=fmodes[[Nt/2+1]]+fmodes[[3*Nt/2+1]];
	(*negative frequencies*)
	Do[faux[[k+1]]=fmodes[[k+Nt+1]],{k,Nt/2+1,Nt-1}];
	faux
];


halve=Compile[{{fmodes,_Complex,1},{Nt,_Integer}},Module[{faux},
	faux=ConstantArray[0.+I*0.,Nt];
	(*positive frequencies*)
	Do[faux[[k+1]]=fmodes[[k+1]],{k,0,IntegerPart[Nt/2-1]}];
	(*sum high-frequency cosines*)
	faux[[IntegerPart[Nt/2+1]]]=fmodes[[IntegerPart[Nt/2+1]]]+fmodes[[IntegerPart[3*Nt/2+1]]];
	(*negative frequencies*)
	Do[faux[[k+1]]=fmodes[[k+Nt+1]],{k,IntegerPart[Nt/2+1],Nt-1}];
	faux
],CompilationTarget->"C"
];


End[];


EndPackage[];
