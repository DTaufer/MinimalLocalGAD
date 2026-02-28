// --------------------------------------------------------------------------------------------------------------------------------------------
// Auxiliary functions
// --------------------------------------------------------------------------------------------------------------------------------------------


function MonomialFactorial( _m )
// Given a monomial m = x^a, it returns a!
	return &*[Factorial(e) : e in Exponents( _m )];
end function;


function Dp( _F )
// Given a homogeneous form F, it returns Fdp, i.e., F written with divided powers
	if IsZero(_F) then return _F; end if; // rule out the F = 0 case
	return &+[_t*MonomialFactorial(_t) : _t in Terms(_F)];
end function;


function Swap( _L, _i, _j )
// Given a list L, and returns the same list after swapping the places of i and j
	local nL, tmp;
	nL := _L;
	tmp := _L[_i];
	nL[_i] := _L[_j];
	nL[_j] := tmp;
	return nL;
end function;


function ListShuffle( _L )
// Given a list L, it shuffles it (= returns a random permutation of L)
	local Ls,n;
	Ls := _L; n := #Ls;
	for i in [0..n-1] do
		Ls := Swap(Ls,Random([1..n-i]),n-i);
	end for;
	return Ls;
end function;


function DualGenerator( _F, _L )
// Given a homogeneous form F and a non-zero linear form L, it returns f_L, i.e. (F_dp)|_{L = 1}.
// ASSUMPTION: L = x0 + a1*x1 + ... + an*xn
	local R,n,K,vars,Fc, names,Raff,dehom;
	R := Parent(_F); n := Rank(R); K := BaseRing(R);
	vars := [ _m : _m in MonomialsOfDegree(R,1)];
	
	Fc := Dp( Evaluate( _F, [2*R.1-_L] cat vars[2..n]) );
	
	names := Names(R);
	Raff := PolynomialRing(K, n-1);
	AssignNames(~Raff, names[2..n]);
	dehom := hom< R -> Raff | [1] cat [ Raff.i : i in [1..(n-1)] ] >;
	
	return dehom(Fc);
end function;


function InvSysMatrix( _fL )
// Given a dual generator fl, it returns the inverse it generates (i.e., the matrix representing all its derivatives)
	local fL,R,n,d,mons;
	
	R := Parent(_fL); n := Rank(R); d := Degree(_fL);
	mons := &join( [ MonomialsOfDegree(R,_d) : _d in [0..d] ]);
	
	return Matrix( BaseRing(R), #mons, [MonomialCoefficient(_fL,_n*_m) : _n,_m in mons] );
end function;


function LocalMeasure( _F, _L )
// Assumption: _L has non-zero x-coeff
	return Rank( InvSysMatrix( DualGenerator( _F, _L ) ) );
end function;


// --------------------------------------------------------------------------------------------------------------------------------------------
// Functions for (B)
// --------------------------------------------------------------------------------------------------------------------------------------------


function MaxEltsFromDeg( _n, _d )
// It returns the vector [N0, ..., Nd], where Ne = maximal number of possibly independent deg-e (equiv. d-e) elements (lines) from the derivatives of an affine degree-d polynomial in n variables
	return [ Min( Binomial(_n+_e-1,_e), Binomial(_n+_d-_e-1,_d-_e) ) : _e in [0.._d] ];
end function;


function BlkSizes( _n, _d )
// It returns a list [N0, ..., Nd], where Ne = the number of deg <= e monomials in n variables
	return [ Binomial(_n+_e,_e) : _e in [0.._d] ];
end function;


function RndBlkDiv( _L, _M, _r )
// Input: block sizes L = [l0, ..., ld], max draws every block M = [m0, ..., md], number of draws r
// Output: a random [r0, ..., rd] such that \sum ri = r and l(i-1) <= ri <= li [convention: l(-1)=0].
// If not possible, it returns [];
	d1 := #_L; // this is d+1
	// CHECK input
	if #_M ne d1 or _L[d1] lt _r or &+_M lt _r then
		"Invalid parameters for RndBlkDiv:";
		#_L, _L;
		#_M, _M;
		_r;
		return [];
	end if;
	
	L0 := [0] cat _L;
	
	// sample r columns - probability density respects the number of monomials of corresponding degree 
	left := _M;
	out_cols := [];
	for i in [1.._r] do
		samplingspace := &cat[ [(L0[_e]+1)..L0[_e+1]] : _e in [1..d1] | not IsZero(left[_e]) ];
		samplingspace := [ _x : _x in samplingspace | _x notin out_cols ];
		//samplingspace;
		draw := Random(samplingspace);
		
		//"drawn", draw;
		drawblock := 1;
		while draw gt _L[drawblock] do
			drawblock +:= 1;
		end while;
		//"drawblock", drawblock;
		left[drawblock] -:= 1;
		//"left", left;
		out_cols cat:= [draw];
	end for;
	
	// sample r rows - drawn randomly from complementary blocks of columns, to avoid trivial zero det
	left := Reverse( [_M[i]-left[i] : i in [1..d1] ] );
	out_rows := {};
	for i in [1..d1] do
		out_rows join:= RandomSubset({(L0[i]+1)..L0[i+1]}, left[i]);
	end for;
	return Sort(SetToSequence(out_rows)), Sort(out_cols);
end function;


// --------------------------------------------------------------------------------------------------------------------------------------------
// Functions for (C)
// --------------------------------------------------------------------------------------------------------------------------------------------


function ContractionLine( _m )
// for a given monomial m = x^{[a,b,c]}, it returns a line of (a+b+c) contraction that brings m to 1
// the output is [ [0,0,0], [1,0,0], [1,1,0], ... [a,b,c] ]
// to get the corresponding monomials, just do Monomial(R,[a,b,c]); 
	local R,n,exps,Lnrd,last,out;
	R := Parent(_m); n := Rank(R);
	exps := Exponents(_m);
	
	Lrnd := ListShuffle( &cat[ [i : j in [1..exps[i]]] : i in [1..n] ] );		// produces a random list containing a-times the index 1 (of a), b-times the index 2 (of b), ...
		
	last := [0: i in [1..n]];
	out := [last];
	for l in Lrnd do
		last[l] +:= 1;
		out cat:= [last];
	end for;
	
	return out;
end function;


function GoodMinor( _r, _fL )
// It returns an (rxr)-minor of its symbolic inverse system matrix of fL.
// It follows contraction chains (method (C)).
	R := Parent(_fL); n := Rank(R); d := Degree(_fL);
	allMons := &cat[ [ _x : _x in MonomialsOfDegree(R,_d) ] : _d in [0..d] ];
	missing := _r;			// counter
	MaxMons := {_m : _m in Monomials(_fL) | Degree(_m) eq Degree(_fL)};
	
	cols := []; 			// the monomials indexing columns
	rows := []; 			// the monomials indexing rows
	
	while missing gt 0 do
		m := Random(MaxMons);
		cl := ContractionLine(m);
		contractions := [ Monomial(R,_l) : _l in cl ];
		monomials := [ m div Monomial(R,_l) : _l in cl ];
		
		for j in [1..#cl] do
			if not (monomials[j] in cols or contractions[j] in rows) then
				cols cat:= [monomials[j]]; rows cat:= [contractions[j]];
				missing -:= 1;
				if missing eq 0 then break j; end if;
			end if;
		end for;
	end while;
	
	return [Index(allMons,_r) : _r in rows], [Index(allMons,_r) : _r in cols];
end function;


// --------------------------------------------------------------------------------------------------------------------------------------------
// MinimalSupports functions (A,B, and C)
// --------------------------------------------------------------------------------------------------------------------------------------------


function MinimalSupportsA( _F, _r )
// Given an integer r and a homogeneous form F, it returns the supports L=x0+a1*x1+...+an*xn such that the natural apolar scheme of F at L has length r
// ASSUMPTION: the number of such points is finite
// Uses method (A): random minors
	R := Parent(_F); n := Rank(R)-1; K := BaseRing(R);
	
	Rsymb<[c]> := PolynomialRing(K,n);
	Rs<[x]> := PolynomialRing(Rsymb,n+1);
	L := x[1] + &+[c[i]*x[i+1] : i in [1..n]];
	names := Names(R);
	AssignNames(~Rs, names);
	
	fL := DualGenerator( Rs!_F, L );
	
	M := InvSysMatrix( fL );
	
	// check that you are not asking for the generic local size
	if Rank(M) le _r then
		"Any generic linear form is minimal, of rank", Rank(M);
		return [];
	end if;	
	
	//"collecting relations";
	Nr := Nrows(M);
	Nc := Ncols(M);
	
	relations := [];
	for _tmp in [1..2*n] do		// a bit of overdetermination
		rows := Sort(Setseq( RandomSubset({1..Nr}, _r+1) ));
		cols := Sort(Setseq( RandomSubset({1..Nc}, _r+1) ));
		relations cat:= [ Minor(M, rows,cols) ];
	end for;
	
	//"constructing the ideal";
	I := ideal<Rsymb | relations >;
	while Dimension(I) gt 0 do
		//Dimension(I);
		rows := Sort(Setseq( RandomSubset({1..Nr}, _r+1) ));
		cols := Sort(Setseq( RandomSubset({1..Nc}, _r+1) )); 
		I +:= ideal<Rsymb | Minor(M, rows,cols) >;	
	end while;
	
	//"Computing points";
	V := Variety(I);
	
	// Now tests every candidate (not all of them are guaranteed to be solutions: they can solve only the chosen minors)
	Ls_sol := {};
	for P in V do
		L := R.1 + &+[P[i-1]*R.i : i in [2..(n+1)]];
		if LocalMeasure( _F, L ) eq _r then 
			Ls_sol join:= { L };
		end if;
	end for;
	return Ls_sol;
end function;


function MinimalSupportsB( _F, _r )
// Given an integer r and a homogeneous form F, it returns the supports L=x0+a1*x1+...+an*xn such that the natural apolar scheme of F at L has length r
// ASSUMPTION: the number of such points is finite
// Uses method (B): block-diagonal minors
	R := Parent(_F); n := Rank(R)-1; K := BaseRing(R); d:= Degree(_F);
	
	Rsymb<[c]> := PolynomialRing(K,n);
	Rs<[x]> := PolynomialRing(Rsymb,n+1);
	L := x[1] + &+[c[i]*x[i+1] : i in [1..n]];
	names := Names(R);
	AssignNames(~Rs, names);
	
	fL := DualGenerator( Rs!_F, L );
	
	M := InvSysMatrix( fL );
	
	// check that you are not asking for the generic local size;
	if Rank(M) le _r then
		"Any generic linear form is minimal, of rank", Rank(M);
		return [];
	end if;
	
	Ls := BlkSizes( n,d );
	Ms := MaxEltsFromDeg( n,d );
	
	relations := [];
	//"collecting relations";
	for _tmp in [1..2*n] do
		rows, cols := RndBlkDiv( Ls, Ms, _r+1 );
		relations cat:= [ Minor(M, rows,cols) ];
	end for;
	
	//"constructing the ideal";
	I := ideal<Rsymb | relations >;
	while Dimension(I) gt 0 do
		//Dimension(I);
		rows, cols := RndBlkDiv( Ls, Ms, _r+1 );
		I +:= ideal<Rsymb | Minor(M, rows,cols) >;	
	end while;
	
	//"Computing points";
	V := Variety(I);
	
	// Now tests every candidate (not all of them are guaranteed to be solutions: they can solve only the chosen minors)
	Ls_sol := {};
	for P in V do
		L := R.1 + &+[P[i-1]*R.i : i in [2..(n+1)]];
		if LocalMeasure( _F, L ) eq _r then 
			Ls_sol join:= { L };
		end if;
	end for;
	return Ls_sol;
end function;


function MinimalSupportsC( _F, _r )
// Given an integer r and a homogeneous form F, it returns the supports L=x0+a1*x1+...+an*xn such that the natural apolar scheme of F at L has length r
// ASSUMPTION: the number of such points is finite
// Uses method (C): contraction chains
	R := Parent(_F); n := Rank(R)-1; K := BaseRing(R);
	
	Rsymb<[c]> := PolynomialRing(K,n);
	Rs<[x]> := PolynomialRing(Rsymb,n+1);
	L := x[1] + &+[c[i]*x[i+1] : i in [1..n]];
	names := Names(R);
	AssignNames(~Rs, names);
	
	fL := DualGenerator( Rs!_F, L );
	
	M := InvSysMatrix( fL );
	
	// check that you are not asking for the generic local size;
	if Rank(M) le _r then
		"Any generic linear form is minimal, of rank", Rank(M);
		return [];
	end if;
	
	//"collecting relations";
	relations := [];
	for _tmp in [1..2*n] do		// a bit of overdetermination
		rows,cols := GoodMinor( _r+1, fL );
		relations cat:= [ Minor(M, rows,cols) ];
	end for;
	
	//"constructing the ideal";
	I := ideal<Rsymb | relations >;
	while Dimension(I) gt 0 do
		//Dimension(I);
		rows, cols := GoodMinor( _r+1, fL );
		I +:= ideal<Rsymb | Minor(M, rows,cols) >;	
	end while;
	
	//"Computing points";
	V := Variety(I);
	
	// Now tests every candidate (not all of them are guaranteed to be solutions: they can solve only the chosen minors)
	Ls_sol := {};
	for P in V do
		L := R.1 + &+[P[i-1]*R.i : i in [2..(n+1)]];
		if LocalMeasure( _F, L ) eq _r then 
			Ls_sol join:= { L };
		end if;
	end for;
	return Ls_sol;
end function;
