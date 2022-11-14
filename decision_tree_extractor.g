LoadPackage("AtlasRep");
LoadPackage("CTblLib");
Read("sr.g");



test1 := function(n)
    local subs, results, g1, g2, sorting_func;

    subs:=AllSmallGroups([2..n]);
    results:=[];
    for g1 in subs do
        for g2 in subs do
            Append(results, [WreathProduct(g1, g2)]);
            # Append(results, [Size(WreathProduct(g1, g2))]);
        od;
    od;
    sorting_func:=function(x, y) return(Size(x)<Size(y));end;;
    Sort(results,sorting_func);
    return(results);
end;


get_entry_props := function(A, B, global, spec, goal_product)
    # Extract props of A, B, goal_product(A, B),
    #    check is_SR(goal_product(A, B)), extend global.
    #
    # Inputs:
    # - ``global`` -- list of props as Union of all entries
    # will be extended here with A, B KnownPropertiesOfObject.
    #
    # - ``spec`` -- list of attributes (funcs) 
    #  to apply (besides KnownPropertiesOfObject) to A, B and goal_product(A, B)
    # 
    # - ``goal_product`` -- funct as WreathProduct, DirectProduct
    #
    # Return:
    # [[A true props], [B true props], [G true props], 
    #  [A false props], [B false props], [G false props], 
    #  [A spec props], [B spec props], [G spec props],
    #  global,
    #  is_SR(goal_product(A, B))]
    # where G is goal_product(A, B)
    #  
    # Example:
    # s3:=Group((1,2), (1,2,3));
    # a4:=AlternatingGroup(4);
    # goal_product := function(A, B)return(WreathProduct(B, A));end;;
    # get_entry_props(s3, a4, [Size, elm->List(ConjugacyClasses(elm, Size))],goal_product);

    local KnownFalsePropertiesOfObject, G, result, props_to_except;

    KnownFalsePropertiesOfObject:=function(A)return(
            Difference(KnownPropertiesOfObject(A),
                       KnownTruePropertiesOfObject(A)));end;;

    props_to_except:=["CanEasilyCompareElements", "CanEasilySortElements","IsAssociative",
                      "IsEmpty", "IsFinite", "IsNonTrivial", "KnowsHowToDecompose",
                      "IsTrivial", "IsSubsetLocallyFiniteGroup"];

    G:=goal_product(A, B);

    result := [];
    Append(result, List([A, B, G], group->Difference(
            KnownTruePropertiesOfObject(group), props_to_except)));
    Append(result, List([A, B, G], group->Difference(
            KnownFalsePropertiesOfObject(group), props_to_except)));
    
    # apply spec props for each of group:
    Append(result, List([A, B, G], group->List(spec, elm->elm(group))));

    # extend global:
    global := Union(global, Union(List([A, B, G],
                group->Difference(KnownPropertiesOfObject(group),
                        props_to_except)))) ;
    Append(result, [global]);
    
    Append(result, [is_SR(G)]);
    return(result);
end;


compute_char:=function(G)
    # compute characters of all direct products DixDj
    # for each irr repr Di, Dj of group G.

    local tbl, xs, len, direct_repr, i, s, j, f;
    tbl:=CharacterTable(G);
    xs := List(Irr(tbl),elm->elm!.ValuesOfClassFunction);
    len := Length(xs);
    direct_repr:=[];
    f :=function(x,y)return(x*y);end;;

    for i in [1..len] do
       # shift for AxB=BxA in indexes;
       s:= i;
       for j in [s..len] do
         # the inner product x[1] * y[1] , x[2] * y[2], ⋯ , x[n] * y[n]:
         Append(direct_repr, [ListN(xs[i], xs[j], f)]);
       od;
    od;
    return(direct_repr);
end;
