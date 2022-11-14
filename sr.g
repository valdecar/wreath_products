
cent:= function(G, g)
        return(Length(List(Centralizer(G, g))));
end;


theta := function(G, g)
     return(Length(Filtered(G, elm1->elm1*elm1=g)));
end;


sum_of_squere:=function(G)
     return(Sum(List(G, elm->theta(G, elm)^3)));
end;


sum_of_cent:=function(G)
     # gap> sum_of_cent(s4);
     # 1032
     return(Sum(List(G, elm->Length(List(Centralizer(G, elm)))^2)));
end;

sum_of_cent1:=function(G)
     # gap> sum_of_cent1(s4);
     # 1032
     return(Sum(List(ConjugacyClasses(G),
                class->Length(List(Centralizer(G,First(class))))^2
                *Length(List(class)))));
end;


is_SR2:=function(G)
        # time = 38104 for $G=$D_{26}\wr S_{2}$

        # complexity:
        # |G|*(2*|G| + 1) = 2*|G|^2+|G| group multiplications
        # |G|+2*|G|= 3*|G| digital multiplications
        

        local squares, cents, thetas_rec, h, square_h, left, right;
       squares:=[];
       cents:=[];
       thetas_rec := rec();
       for h in G do
          # FOR theta:
          # some elements g has no solution of h^2=g
          # so We force them to thetas_rec with zero value:
          if not IsBound(thetas_rec.(String(h))) then
            thetas_rec.(String(h)) := 0;
          fi;

          # for others collect square solutions counts:
          square_h:=h^2;  # 1 group multiplication
          if not IsBound(thetas_rec.(String(square_h))) then
             thetas_rec.(String(square_h)) := 1;
          else
             thetas_rec.(String(square_h)) := thetas_rec.(String(square_h))+1;
          fi;
          # END FOR

          # Add(squares, h^2);
          # 2*|G| group multiplications
          # 1 digital multiplications 
          Add(cents, Size(Centralizer(G, h))^2);
        od;
        left:=Sum(cents);
        # |G|*2 digital multiplications:
        right:=Sum(G, e->thetas_rec.(String(e))^3);
        Print("left: ", left, "\n");
        Print("right: ", right, "\n");
        return(left=right);
        # return(Sum(cents)=Sum(Collected(squares), e->e[2]^3));
end;


is_SR1:=function(G)
    # gap> is_SR1(WreathProduct(ElementaryAbelianGroup(2^4),s2));

    # gap> w:=WreathProduct(s3,s2);
    # gap> is_SR(w);
    # true
    local left, right;
    left:= sum_of_squere(G);
    Print("left: ", left, "\n");
    right:=sum_of_cent(G);
    Print("right: ", right, "\n");
    return(left=right);
end;


compute_reduction := function(G)
    local tbl, xs, len, direct_repr, i, s, j, k, f, cls_sizes, tmp, tmp1;
    tbl:=CharacterTable(G);
    xs := List(Irr(tbl),elm->elm!.ValuesOfClassFunction);
    cls_sizes := SizesConjugacyClasses(tbl);

    len := Length(xs);
    Print("Length xs ~ count of characters of G = ",len," done","\n");
    # Print(Concatenation("Length xs ~ count of characters of G = ",len," done"),"\n");

    direct_repr:=[];
    f :=function(x,y,z,t)return(x*y*z*t);end;;

    for k in [1..len] do
     tmp:=[];
     for i in [1..len] do
      # shift for AxB=BxA in indexes;
      s:= i;
      for j in [s..len] do
       # the inner product (class_size[i] * xs[i])*xs[j]*xs[k]:
       tmp1:=Sum(ListN(cls_sizes, xs[k], xs[i], xs[j], f))/Size(G);
       Append(direct_repr, [tmp1]);
       if tmp1 > 1 then
        Print("i,j,k = ",i," ",j," ", k, "\n");
        Print("\n count of such  = ", tmp1, "\n");
        Print("xs[k] = ", xs[k], "\n");
        Print("xs[i] = ", xs[i], "\n");
        Print("xs[j] = ", xs[j], "\n");
        # Print("sum = ", tmp1, " xs[k][e] = ", xs[k][1], " xs[j][e] = ", xs[j][1], " xs[i][e] = ", xs[i][1], "\n");
         
       fi;
       # in that case result for g1:=SmallGroup(6,1); will be:
       # [ [ 1, 0, 0, 1, 0, 1 ], [ 0, 1, 0, 0, 0, 1 ], [ 0, 0, 1, 0, 1, 1 ] ] 
       # in p. 200 of Hamermesh [A1,A2,E] = [direct_repr[1][-1], direct_repr[2][-1], direct_repr[3][-1]]; 
       # Append(tmp, [Sum(ListN(cls_sizes, xs[k], xs[i], xs[j], f))/Size(G)]);         
        
       od;
    od;
    # Append(direct_repr, [tmp]);
    
    Print("progress: k = ",k," from ", len,"\n");
    # Print(Concatenation("progress: k = ",k," from ", len),"\n");
    
    od;
    return(direct_repr);
end;

is_SR:=function(G)
        # time = 23376 for G=$D_{26}\wr S_{2}$
        # check if group G is Simply Reducable
        return(Length(Filtered(compute_reduction(G), elm->elm>1))=0);
end;

