test_S4xA5_s2:=function(S4xA5, subs2)
    
    # gap> subs2:=Filtered(subsc, elm->Size(elm[1])=2);
    # 
    # this is long:
    # gap> subsc:=ConjugacyClassesSubgroups(au);
    
    # or use
    # cls_au_flat:=Flat(List(cls_au_filt, cl->List(cl)));
    
    # gap> au:=AutomorphismGroup(S4xA5);
    
    local s2, homs, sems, sem, res, is_sr, irrs;
    s2:=Group((1,2));

    # homs:=List(subs2, sub-> GroupHomomorphismByImages(s2, sub, [(1,2)], GeneratorsOfGroup(sub)));
    homs:=List(subs2, sub-> GroupHomomorphismByImages(s2, sub[1], [(1,2)], GeneratorsOfGroup(sub[1])));
    sems:=List(homs, hom->SemidirectProduct(s2, hom, S4xA5));
    res:=[];
    irrs:=[];
    for sem in sems do
        #is_sr:=is_SR(sem);
        #if is_sr then
        #    Add(res, [sem, is_sr]);
            Add(irrs, [sem, Collected(List(Irr(sem), e->e[1]))]);
        #fi;
    od;
    return(irrs);
    # return(res);
end;
