testWreath := function(G, H)
    local
        shiftL, shiftR,
       
	GP,	    # preimage of <G> (if homomorphism)
	Ggens,Igens,#permutation images of generators
    	alpha,      # action homomorphism for <H>
	permimpr,   # product is pure permutation groups imprimitive
        I,          # image of <alpha>
        grp,        # wreath product of <G> and <H>, result
        gens,       # generators of the wreath product
	gens1,      # generators of first base part
        gen,        # one generator
        domG,       # domain of operation of <G>
        degG,       # degree of <G>
        domI,       # domain of operation of <I>
        degI,       # degree of <I>
        shift,      # permutation permuting the blocks
	perms,      # component permutating permutations
	basegens,   # generators of base subgroup
	hgens,	# complement generators
	components, # components (points) of base group
        rans,       # list of arguments that have '.sCO.random'
	info,	# info record
        i, k, l;    # loop variables

    domG := MovedPoints( G );
    Ggens:=GeneratorsOfGroup(G);
    degG := Length( domG );
    
    domI := MovedPoints( H );
    Igens:=GeneratorsOfGroup(H);
    domI := MakeImmutable(Set(domI));
    degI := Length( domI );
    hgens:=[];
    gens := [];
    for gen  in Igens  do
        Print("gen=",gen, "\n");
        shift := [];
	shiftL:=[];
	shiftR:=[];
        for i  in [1..degI]  do
            k := Position( domI, domI[i]^gen );
            Print("i=", i, " k=", k, "\n");
	    for l  in [1..degG]  do
	        Print(" l=", l, "\n");
                Print("(i-1)*degG+l = ", (i-1)*degG+l, "\n");
		Print("(k-1)*degG+l = ", (k-1)*degG+l, "\n");
		Add(shiftL, (i-1)*degG+l);
		Add(shiftR, (k-1)*degG+l);
		shift[(i-1)*degG+l] := (k-1)*degG+l;
            od;
        od;
        Print("gen=",gen, "\n");
	Print("shiftL:", shiftL, "\n");
	Print("shiftR:", shiftR, "\n");
        Print("shift:", shift);
	Print("\n");
	shift:=PermList(shift);
        Add( gens, shift );
	Add(hgens, shift );
    od;
    return(hgens);
end;
