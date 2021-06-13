#!/usr/bin/env wolframscript

(*The following are more important modules in the Mathematica language 
in which we programmed the calculation of the critical values of the parameter t.*)

cellsFromLeft::usage="cellsFromLeft[l0_,r0_,h0_,q0_] 
computes the minimum number of cells in a layer, 
such that the last one is suboptimal.";
cellsFromLeft[l0_,r0_,h0_,q0_]:=Module[
	{lf,l=l0,r=r0,h=h0,q=q0, f1, f2,alist},
	f1[x_]:=q*x; f2[x_]:=h(*@-@*)(h(*@-@*)x)/q;
	alist={l};(*there will always be at least one cell*)
	While[(l=f1[l])<h/2,AppendTo[alist,l]];
	lf=Last[alist];(*boundary just before the middle*)
	l=Min[f1[lf],f2[lf]];(*boundary just after the middle*)
	If[l<r0,AppendTo[alist,l]];
	(*otherwise it will only be a single cell*)
	While[(l=f2[l])<r0,AppendTo[alist,l]];
	Length[alist]]

findSubcritTs::usage="findSubcritTs[t00_,t01_, fj_, l0_,h0_,q0_] 
returns a list of subcritical t's, that is where the definition 
of L(t) is changing";
findSubcritTs[t00_,t01_, fj_, l0_,h0_,q0_]:=Module[
	{h = h0, q = q0,numOfSubcritTs,max,min,t0=t00,t1=t01, l0Inv},
	max=fj[t0];min=fj[t1];
	numOfSubcritTs=max(*@-@*)min;
	If[numOfSubcritTs>0,
		Sort[Union[
		l0Inv[y_]:=InverseFunction[l0][y];
		Map[l0Inv[h/(q^(max(*@-@*)#)*(q+1))]&,Range[0,
		numOfSubcritTs(*@-@*)1]],
		{t0,t1}]],
		{t0,t1}]]

findCritTsFromFunctions::usage=
"findCritTsFromFunctions[subIntBLX_, h0_] process the layers 
on the observed subinterval between a pair of adjacent critical
t's sct0 and sct1 where fj has jump discontinuities. There we 
know that fj is constant fj(t) = fj((sct0+sct1)/2) =: f0, for 
every t on the subinterval (sct0, sct1).  We can also compute 
fk[sct0] =: k0. Thus, we have L_{k(*@-@*)1} uniquely defined. 
We are interested in whether there is a critical t on the 
observed subinterval where fk has jump discontinuities. 
The rest of the story is kept in the comments.";
findCritTsFromFunctions[subIntBLX_, h0_]:=Module[
	{subIntBL = subIntBLX, h = h0, q, f1, f2, lkMinOne, 
	fk, fj, t0, t1, p,max, min,numOfSubcritTs,subcritTs,
	L0Inv, i,j, l0, r0, k0, j0,sct0,sct1,critTsList,t,ct},
	critTsList={};
	h = h0;
	q = (1+3(h+1/3))/(3(h+1/3));
	f1[x_]:=q*x; f2[x_]:=h(*@-@*)(h(*@-@*)x)/q;
	lkMinOne[t_,k0_, j0_]:= Nest[f2,Nest[f1,l0[t],j0+1], k0(*@-@*)j0(*@-@*)2]; 
	For[i=1,i<=Length[subIntBL],i++,
		p=Part[subIntBL,i];
		t0=Part[p,1];t1=Part[p,2];
		l0[t_]:=Part[p,3][t,h];
		r0[t_]:=Part[p,4][t,h];
		fk[t_]:=cellsFromLeft[l0[t],r0[t],h,q];
		fj[t_]:=Floor[Log[h/(l0[t]*(q+1))]/Log[q]];
		subcritTs = findSubcritTs[t0,t1, fj, l0,h,q];
		For[j=1,j<Length[subcritTs],j++,
(*Critical t exists and fk(t) decreases for one, 
if L_{k(*@-@*)1} reaches R_{0} or critical t exists and fk(t) 
increases for one, if R_{0} surpasses L_ {k}. 
In the second case there are two possibilities, 
if k(*@-@*)1 > j0 we are using second definition for L_{k} 
and have to check if R_{0} > f_{2}(L_{k(*@-@*)1}), 
otherwise we have to check if 
R_{0} > f_{1}(L_{k(*@-@*)1}).*)
			sct0=Part[subcritTs,j];
			sct1=Part[subcritTs,j+1];
			k0=fk[sct0];
			j0=fj[(sct0+sct1)/2];
(*L_{0}, f_{1}, f_{2} and R_{0} are linear functions. L_{k} is the
composition of k linear functions, so it is also a linear function 
for every fixed k. For all the described cases we have to find the 
intersection of two lines. First we check whether the two lines are 
parallel. If not, we can check for intersection. Finaly we check if 
this happens on the observed subinterval.*)
			If[(D[lkMinOne[t,k0,j0],t]/.t(*@-@*)>t0)
			!=(D[r0[t],t]/.t(*@-@*)>t0),
				ct=Solve[lkMinOne[t,k0,j0]==r0[t],t];
				ct = t /. ct[[1]];
				If[sct0<ct&&ct<sct1,AppendTo[critTsList,ct]]
			](*Derivative[lkM] > 0?*)
			If[k0(*@-@*)1>j0,
				If[(D[f2[lkMinOne[t,k0,j0]],t]/.t(*@-@*)>t0)!=
					(D[r0[t],t]/.t(*@-@*)>t0),
					ct=Solve[f2[lkMinOne[t,k0,j0]]
					==r0[t],t];
					ct = t /. ct[[1]];
					If[sct0<ct&&ct<sct1,
					AppendTo[critTsList,ct]]],
				If[(D[f1[lkMinOne[t,k0,j0]],t]/.t(*@-@*)>t0)
				!=(D[r0[t],t]/.t(*@-@*)>t0),
					ct=Solve[f1[lkMinOne[t,k0,j0]]
					==r0[t],t];
					ct = t /. ct[[1]];
					If[sct0<ct&&ct<sct1,
					AppendTo[critTsList,ct]]
				];
			];
		];
(*If on the observed subinterval L_{k(*@-@*)1} does not intersect R_{0} 
and also R_{0} does not intersect L_{k}, than that means that 
there are no critical t's on the observed subinterval.*);
	];
critTsList]


findCritTsBottomLayers[]:=Module[
	{bndStubL,bndStubR,bndFrameL,bndFrameR,
	bndRayL,bndRayR,subIntBL2,subIntBL3,t,h},
	bndStubL[t_,h_] := (h(1+((*@-@*)2+3 h) t))/(2+3 h);
	bndStubR[t_,h_] := h (*@-@*) bndStubL[1 (*@-@*) t,h];
	bndFrameL[t_,h_] := (h(*@-@*)1)*t;
	bndFrameR[t_,h_] := h (*@-@*) bndFrameL[1(*@-@*)t,h];
	bndRayL[t_,h_]:= h/4;
	bndRayR[t_, h_] := (3*h)/4;
	subIntBL2 := { 
	{0, 1/4, bndStubL, bndFrameR}, 
	{1/4, 3/8, bndRayL, bndFrameR}, 
	{3/8, 1/2, bndRayL, bndRayR}};
	subIntBL3 := {
	{0, 1/4, bndStubL, bndFrameR},
	{1/4, 3/8, bndRayL, bndFrameR},
	{3/8, 1/2, bndFrameL, bndFrameR}};
{findCritTsFromFunctions[subIntBL2, 4/3],
findCritTsFromFunctions[subIntBL3, 5/3]}]


findSftIntFromFunctions::usage=
"findSftIntFromFunctions[subIntTLX_, h0_,eps0_] estimate optimal
critical t's";
findSftIntFromFunctions[subIntTLX_, h0_,eps0_]:=Module[
	{subIntTL = subIntTLX, h = h0,t0, t1,l0,
	r0,p,fq,fk,sftIntList,i,eps=eps0},
	sftIntList={};
	For[i=1,i<=Length[subIntTL],i++,
		p=Part[subIntTL,i];
		t0=Part[p,1];t1=Part[p,2];
		l0[t_]:=Part[p,3][t];
		r0[t_]:=Part[p,4][t];
		fq[t_] := Part[p,5][t];
		fk[t_]:=cellsFromLeft[l0[t],r0[t],h, fq[t]];
		AppendTo[sftIntList,findSftIntBisection[t0,t1,fk,eps]];];
sftIntList]


findSftIntBisection[t00_, t01_, fk_, eps0_]:= Module[
	{t0 = t00, t1 = t01, min, max,
	numOfCTs,i,a,b,c, eps = eps0,sftSubintList},
	sftSubintList = {};
	min=fk[t0];
	max=fk[t1];
	numOfCTs=max(*@-@*)min;
	For[i=1,i<=numOfCTs,i++,
		a=t0;b=t1;
		c=(a+b)/2;
		While[b(*@-@*)a>eps,
			If[fk[c](*@-@*)(max(*@-@*)i)<=0,a=c,b=c];
			c=(a+b)/2]; (*end while*)
		AppendTo[sftSubintList,{a,b}]];(*end for i*)
	Sort[sftSubintList]]


findSftIntTopLayers[]:=Module[
	{eps,bndRayL1,bndRayR1,fq1,bndFrameL2,bndFrameR2,fq2i,fq2ii,
	bndFrameL3,bndFrameR3,fq3,subIntTL1,subIntTL2,subIntTL3},
	eps = 10^(*@-@*)10;
	(*TL1*)
	bndRayL1[t_]:= 1/4; 
	bndRayR1[t_]:= 3/4;
	fq1[t_] := 4/3(*@-@*)t/4;
	(*TL2*)
	bndFrameL2[t_] := 1/(3 t);
	bndFrameR2[t_]:= 1;
	fq2i[t_] := 1+1/(3+3 t);(*first part*)
	fq2ii[t_] :=  6/5; (*second part*)
	(*TL3*)
	bndFrameL3[t_] := 2/(3 t);
	bndFrameR3[t_] := 1;
	fq3[t_] := 1+1/(3+3 t);
	(*{t0, t1, l0, r0, fq}*)
	subIntTL1 := { 
	{0, 1/3, bndRayL1,bndRayR1,fq1}};
	subIntTL2:= {
	{1/3, 2/3,bndFrameL2, bndFrameR2,fq2i},
	{2/3, 1, bndFrameL2,bndFrameR2,fq2ii}};
	subIntTL3 := {
	{2/3, 1, bndFrameL3, bndFrameR3, fq3}};
{findSftIntFromFunctions[subIntTL1,1,eps],
findSftIntFromFunctions[subIntTL2, 4/3,eps],
findSftIntFromFunctions[subIntTL3, 5/3,eps]}]


disjointIntervalsQ[intervals0_]:= Module[
{intervals = intervals0,numOfIntervals,aibi,order, result},
	numOfIntervals = Length[intervals];
	intervals=Flatten[intervals];
	aibi = Table[{i,i}, {i, 1, numOfIntervals}];
	aibi = Flatten[aibi];
	order = Ordering[intervals]; (*tells to which place each one 
	has to go*)
	result = aibi[[order]]; (*we place indexes in this order*)
	result = makePairs[result]; (*back to intervals*)
	Total[Map[If[Part[#,1] == Part[#,2], 1, 0]&, result]]
	 == numOfIntervals (*if they are the same number in each pair*) 
	&& Length[DeleteDuplicates[intervals]] 
	== Length[intervals] (*intervals must not touch*)
]


makePairs[checklist0_]:=
Module[{checklist = checklist0,checklist2, ai, bi},
	checklist2 = {};
	For[i = 1, i < Length[checklist], i = i + 2,
	ai = Part[checklist, i];
	bi = Part[checklist, i+1];
	AppendTo[checklist2,{ai,bi}]
];
checklist2]


inflateTs[tsList_,eps_]:=Module[
	{inflatedTsList,ais,bis},
	ais =tsList+eps; (*we blow up exact critical t's*)
	bis=tsList(*@-@*)eps; (*for eps*)
	Inner[List, ais, bis, List]]


disjointSafetyIntervalsQ[]:=Module[
{sftIntBL,sftIntBLSim,sftIntTL,
sftIntTLSim,sftIntBH,sftIntBHSim,sftIntAll},
	sftIntBL=inflateTs[Flatten[findCritTsBottomLayers[]], 10^(*@-@*)10];
	sftIntBLSim = 1(*@-@*)sftIntBL;
	sftIntTL = Flatten[findSftIntTopLayers[],2];
	sftIntTLSim = 1(*@-@*)sftIntTL;
	sftIntBH = findSftIntBlakHoles[];
	sftIntBHSim = 1(*@-@*)sftIntBH;
	sftIntAll = Union[sftIntBL,sftIntBLSim,
	sftIntTL,sftIntTLSim,sftIntBH,sftIntBHSim];
	disjointIntervalsQ[sftIntAll]]


findPartsOfGlobalFunction[]:=
Module[
	{fklist,tflist,i,tlist,flist},
	fklist={};
	For[i=1,i<Length[tflist],i=i+2,
	tlist =Part[tflist,i];
	flist =Part[tflist,i+1];
	AppendTo[fklist, Piecewise[
	Table[{Part[flist,i],
	Part[tlist,i]<= t<= Part[tlist,i+1]},{i,1,Length[flist]}],0]];
];
fklist]


midRationals[sftIntAll_] := Module[
{midRats,allbis,allais,i},
	midRats = {};
	allbis = Map[Max,sftIntAll];
	allbis = Part[allbis,1;;Length[allbis](*@-@*)1];
	allais = Map[Min,sftIntAll];
	allais = Part[allais,2;;Length[allais]];
	Table[Rationalize[N[Mean[{allbis[[i]], allais[[i]]}]], 
	(allais[[i]](*@-@*) allbis[[i]])/2],{i, 1,Length[allais]}]];


findRationalsForEvaluation[]:=Module[
{subCritTs,sftIntSubCrit,sftIntBL,sftIntBLSim,sftIntTL,sftIntTLSim,
sftIntBH,sftIntBHSim,sftIntAll,allbis,fk,fklist,fkGlobal,midRats},
	subCritTs = Union[{0, 1/4, 3/8, 1/2},1(*@-@*){0, 1/4, 3/8, 1/2},
	{0,1/3,2/3,1},
	{0,3/8,5/8},1(*@-@*){0,3/8,5/8}];
	sftIntSubCrit = inflateTs[subCritTs, 10^(*@-@*)10];
	sftIntBL=inflateTs[Flatten[findCritTsBottomLayers[]], 10^(*@-@*)10]; 
	(*with Map[Mean,sftIntBL] we get exact ct's back*)
	sftIntBLSim = 1(*@-@*)sftIntBL;
	sftIntTL = Flatten[findSftIntTopLayers[], 2];
	sftIntTLSim = 1(*@-@*)sftIntTL;
	sftIntBH = findSftIntBlakHoles[];
	sftIntBHSim = 1(*@-@*)sftIntBH;
	sftIntAll = Sort[Union[sftIntSubCrit,sftIntBL,
	sftIntBLSim,sftIntTL,sftIntTLSim,sftIntBH,sftIntBHSim]];
	midRats = midRationals[sftIntAll]]



(*The following is an example of using modules. On all systems of safety intervals, 
we add the functional values of all functions, add constants and calculate the maximum. 
The result is an estimate for the upper bound on the number of basis points 
of a partial drawing of a complete graph.*)

fklist = findPartsOfGlobalFunction[];
(*(B_2 + B_3 + L_1 + L_2 + L_3 + P_1 + P_2)[t]*)
fk[t_]=fklist[[1]]+fklist[[2]]+fklist[[3]]+
fklist[[4]]+fklist[[5]]+fklist[[6]]+fklist[[7]];
fkGlobal[t_]:=fk[t]+fk[1 (*@-@*) t]+39+7;
rats = findRationalsForEvaluation[]
Max[Map[fkGlobal,rats]]

