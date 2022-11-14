find_represent:=function(G, char)
    # char := ,First(Irr(G),e->e[1]=deg);
    # find_represent(d1,First(Irr(G),e->e[1]=deg));
    local hom, cls;
    hom:=IrreducibleRepresentationsDixon(G,char:unary);
    cls:=List(ConjugacyClasses(G),Representative);
    Print("classes representetives: ", cls, "\n");
    Print("matrix Represetations of char irr representation:\n");
    return(List(cls,e->Image(hom,e)));
end;
