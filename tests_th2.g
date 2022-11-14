Read("sr.g");

idxs1 := H-> Cartesian(H, H);

# (u,v) in HxH: u in v^H:
idxs2 := H->Flat(List(List(H, f->f^H),hs->List(hs)));
idxs21 := H->List(H, f->[f, f^H]);


# (u,v) in HxH: u not in v^H:
idxs3 := H->List(H, f->[f, Filtered(H, h->not (h in f^H))]);

########################

test_l4:=function(H, sigma)

    # gap> test_l4(s3, (1, 4)(2, 5)(3, 6));
        
    local w, l1, l2, l3;
    w:=WreathProduct(H, Group((1,2)));
    l1 := Sum(idxs1(H), uv->theta(w, sigma * uv[1] * uv[2]^sigma));
    Print("statement 1: ", l1, "\n");

    # u, v are not conjugate:
    l2:=ForAll(idxs3(H), vu->ForAll(vu[2], u->theta(w, vu[1]*u^sigma)=theta(H, vu[1])*theta(H, u)));
    Print("statement 2:", l2, "\n");

    # u, v are conjugate:
    l3:=ForAll(idxs21(H), vu->ForAll(vu[2], u->theta(w, vu[1]*u^sigma)=theta(H, vu[1])*theta(H, u)+Size(Centralizer(H, u))));
    Print("statement 3:", l3);
    return((l1=0) and l2 and l3);
end;

########################
test_sum_8_left:=function(H)

        local w;
        w:=WreathProduct(H, Group((1,2)));
        return(Sum(w, g->theta(w, g)^3));
end;

test_sum_8_right:=function(H)
       # gap> test_sum_8_right(s3);
       # right: 11232

       local r1, r2, right, right1;
       r1:=Sum(idxs3(H), fh->theta(H, fh[1])^3 * Sum(fh[2], h->theta(H, h)^3));
       # \sum_{(u,v) \in idxs21}(t(u)*t(v)+|C(u)|)^3:
       r2:= Sum(idxs21(H), vu->Sum(vu[2], u->(theta(H, u)*theta(H, vu[1])+Size(Centralizer(H, u)))^3));

       right:=r1+r2;
       Print("right: ", right, "\n");

       r1:=Sum(idxs1(H), uv->Product(uv, x->theta(H, x)^3));
       r2:=Sum(idxs2(H), u->(3*theta(H, u)^4*Size(Centralizer(H, u))+3*theta(H, u)^2*Size(Centralizer(H, u))^2+Size(Centralizer(H, u))^3));
       right1:=r1+r2;
       Print("right1: ", right1);
       return(right1);
end;

test_sum_5_left:=function(H)
       local w;
       w:=WreathProduct(H, Group((1,2)));
       return(Sum(w, g->Size(Centralizer(w, g))^2));
end;

test_sum_5_right:=function(H)
       # gap> test_sum_5_right(s3);
       # rigth1:11232
       # gap> test_sum_5_left(s3);
       # 11232

       local r1, r2, r3, right1, right2;

       # \sum_{f, h \in H}(4*|C_{H}(h*f)|^2):
       r1 := Sum(idxs1(H), hf->4*Size(Centralizer(H, Product(hf)))^2);

       # sum_{(h,f)\in HxH| h \in f^H}4*|C_{H}(h)|^4:
       r2 := Sum(idxs2(H), h->4*Size(Centralizer(H, h))^4);

       # sum_{(h,f)\in HxH| h \not \in f^H}|C_{H}(h)|^2*|C_{H}(f)|^2:
       r3 := Sum(idxs3(H), fh->Size(Centralizer(H, fh[1]))^2 * Sum(fh[2], h->Size(Centralizer(H, h))^2));

       right1:=r1+r2+r3;
       Print("rigth1:", right1, "\n");

       r1 := 4*Size(H)*Sum(H, h->Size(Centralizer(H, h))^2);
       r2 := Sum(idxs1(H), hf->Product(hf, x->Size(Centralizer(H, x))^2));
       r3 := 3*Sum(idxs2(H), h->Size(Centralizer(H, h))^4);
       right2:=r1+r2+r3;
       Print("rigth2:", right2);
       return(right2);
end;


test_sum_8_1 := function(H)
        # gap> test_sum_8_1(Group((1,2),(1,2,3)));
        # true
        #
        # gap> test_sum_8_1(Group((1,2),(1,2,3,4)));
        # true

        local idxs, left, right;
        idxs:=Cartesian(H, H);
        left:=Sum(idxs, uv->Product(uv, x->theta(H, x)^3));
        Print("left: ", left,"\n");
        right:=Sum(idxs, uv->Product(uv,x->Size(Centralizer(H, x))^2));
        Print("right: ", right, "\n");
        return(left=right);
end;

####################################

l3_left1:=function(H, h, f, sigma)
    # wreath product with S2 only!
    # gap> Length(l3_left2(s4, (1,2,3,4),(1,2,3), (1, 5)(2, 6)(3, 7)(4, 8)));
    # 4

    # local sigma;
    # sigma := (1, 4)(2, 5)(3, 6);
    return(List(Centralizer(H, h*f), y-> y^h * y^sigma));
end;

l3_left2:=function(H, h, f, sigma)
    # wreath product with S2 only!
    # gap> Length(l3_left2(s4, (1,2,3,4),(1,2,3), (1, 5)(2, 6)(3, 7)(4, 8)));
    # 4

    # local sigma;
    # sigma := (1, 4)(2, 5)(3, 6);
    return(List(Centralizer(H, f*h), t-> sigma * (f^(-1))*t * (t*h^(-1))^sigma));
end;

l3_right:=function(w, h, f, sigma)
   # wreath product with S2 only!
   # local sigma;
   # sigma := (1, 4)(2, 5)(3, 6);
    
   return(Centralizer(w, sigma * h * f^sigma ));
end;

test_l3:=function(H, sigma)
   # gap> test_l3(s3, (1, 4)(2, 5)(3, 6));
   # true
   # gap> test_l3(s2, (1, 3)(2, 4));
   # true
   # gap> test_l3(s4, (1, 5)(2, 6)(3, 7)(4, 8));
   # true
   # now one can trust lemma 3!

   local h, f, s2, r, r1, l1, l2, l3, test;

   s2 := Group((1,2));
   test := true;

   for h in H do
      for f in H do 
         r:=l3_right(WreathProduct(H, s2), h, f, sigma);
         r1:=Size(r);
         l1:=Length(l3_left1(H, h, f, sigma));
         l2:=Length(l3_left2(H, h, f, sigma));
         l3:=Size(Centralizer(H, h*f));

         test:= test and (r1=l1+l2) and (r1=2*l3);

        if not ((r1=l1+l2) and (r1=2*l3)) then
            Print("r1,l1,l2 :", r1, " ", l1, " ", l2,"\n");
            Print("r1=l1+l2: ", r1=l1+l2,"\n");
            Print("h: ", h,"\n");
            Print("f: ", f,"\n","\n");
        fi;
      od;
   od;
   return(test);
end;

l1_left1:=function(H, h, f, v, sigma)
    # wreath product with S2 only!
    local result, t, u, first, second;
    # sigma := (1, 4)(2, 5)(3, 6);
    result:=[];
    for t in Centralizer(H, h) do
        for u in Centralizer(H, f) do
            first := t * u^sigma;
            if not (first in result) then
                Add(result, first);
            fi;
            second := sigma * u*v * (t*(v^(-1)))^sigma;
            if not (second in result) then
                Add(result, second);
            fi;

        od;
    od;
    return(result);
end;


test_l1 := function(H, sigma)

    # gap> test_l1(s2, (1, 3)(2, 4));
    # true

    # gap> test_l1(s3, (1, 4)(2, 5)(3, 6));
    # true

    # gap> test_l1(s4, (1, 5)(2, 6)(3, 7)(4, 8));
    # true

    # now one can trust lemma 1!

    local s2, f, v, h, r, r1, l1, l2, test;
    test := true;
    s2 := Group((1,2));

    for f in H do
        for v in H do
            h:=f^v;
            r:= Centralizer(WreathProduct(H, s2), h * f^sigma);
            r1:=Size(r);
            l1:=Length(l1_left1(H, h, f, v, sigma));
            l2:=Size(Centralizer(H, h));
            if not ((r1=l1) and (r1=2*l2^2)) then
                Print("r1,l1:", r1, " ", l1, "\n");
                Print("h: ", h,"\n");
                Print("f: ", f,"\n");
                Print("v: ", v,"\n","\n");
            fi;
            
            test := test and (r1=l1) and (r1=2*l2^2);
        od;
    od;
    return(test);
end;


l2_left:=function(H, h, f, sigma)
    # wreath product with S2 only!
    local result, x, y, first;
    result:=[];
    for x in Centralizer(H, h) do
        for y in Centralizer(H, f) do
            first := x * y^sigma;
            if not (first in result) then
                Add(result, first);
            fi;
        od;
    od;
    return(result);
end;


test_l2:=function(H, sigma)
      # gap> test_l2(s2, (1, 3)(2, 4));
      # true

      # gap> test_l2(s3, (1, 4)(2, 5)(3, 6));
      # true

      # gap> test_l2(s4, (1, 5)(2, 6)(3, 7)(4, 8));
      # true

      # now one can trust lemma 2!

      local s2, f, h, r, r1, l1, l2, l3, test;
      test:=true;
      s2:=Group((1,2));
      for f in H do
        for h in H do
            if not (h in f^H) then
                r:= Centralizer(WreathProduct(H, s2), h * f^sigma);
                r1:=Size(r);
                l1:=Length(l2_left(H, h, f, sigma));
                l2:= Size(Centralizer(H, h));
                l3:= Size(Centralizer(H, f));
                if not ((r1=l1) and (r1=l2*l3)) then
                    Print("r1,l1:", r1, " ", l1, "\n");
                    Print("h: ", h,"\n");
                    Print("f: ", f,"\n");
                    Print("l2,l3:", l2, " ", l3, "\n", "\n");
                fi;
                
                test := test and (r1=l1) and (r1=l2*l3);
            fi;
        od;
    od;
    return(test);
end;


############################
find_conjugates1:=function(H)
     # seems correct,
     # used for testing find_conjugates
     return(List(ConjugacyClasses(H),
                    class->Length(List(Centralizer(H, First(class))))
                        *Length(List(class))));
end;


find_conjugates:=function(H)
     # gap> Sum(find_conjugates(s4));
     # 120
     # gap> Sum(find_conjugates(s3));
     # 18
     
     local result, v, u, tmp;
     result:=[];
     tmp:=[];
     for v in H do
        for u in v^H do
          if not (u in tmp) then
            Add(tmp, u);
            Add(result, Length(List(Centralizer(H,u))));
          fi;
        od;
    od;
    # return(result);
    return(tmp);
end;



kolesnikov_right := function(H)
    # gap> kolesnikov_right(ElementaryAbelianGroup(2^3));
    # 114688 # == 2^(4*3)*(4+3*2^3);
    # gap> kolesnikov_right(ElementaryAbelianGroup(2^4));
    # 3407872 # == 2^(4*4)*(4+3*2^4);
    # gap> kolesnikov_right(Group((1,2),(1,2,3)));
    # 6876

    local idxs, right, left;

    # (u,v) in HxH: u in v^H:
    idxs := Flat(List(List(H, f->f^H),hs->List(hs)));

    right := 3*Sum(idxs, h->Size(Centralizer(H, h))^4);
    # right:=3*Sum(List(ConjugacyClasses(H),
    #             class->Length(List(Centralizer(H,First(class))))^4
    #             *Length(List(class))));
    left:= 4*Size(H)*Sum(List(H, elm->cent(H, elm)^2));
    return(left+right);
end;



kolesnikov_left := function(H)
    # gap> kolesnikov_left(ElementaryAbelianGroup(2^3));
    # 114688 # == 2^(4*3)*(4+3*2^3);
    # gap> kolesnikov_left(ElementaryAbelianGroup(2^4));
    # 3407872 # == 2^(4*4)*(4+3*2^4);
    # gap> kolesnikov_left(Group((1,2),(1,2,3)));
    # 6876

    local idxs, thetas;
    # idxs:=find_conjugates(H);

    # (u,v) in HxH: u in v^H:
    idxs := Flat(List(List(H, f->f^H),hs->List(hs)));

    thetas:=List(idxs, elm->[elm, theta(H, elm)]);
    #Print("thetas: ", thetas);
    return(Sum(List(thetas,
            elm->(3*elm[2]^4*cent(H, elm[1])
            + 3*elm[2]^2*cent(H, elm[1])^2
            + cent(H, elm[1])^3))));
end;


kolesnikov_SR_sufficient1:=function(H)

   # time = 228 for $G=$D_{26}\wr S_{2}$

   # complexity:
   # 4*|H|^2+|H| group multiplications
   # 16*|H|^2 + |H| + 2 digital multiplications

   local steps, cents, cents_tmp, theta_counter, thetas_rec, h1, square_h1,
    Hs, h2, cs_h2_H, left, right; 
   steps:=[];
   cents:=[];
   cents_tmp:=rec();
   
   theta_counter:=0;
   thetas_rec:=rec();
   
   # |H|*(2*|H|+2*|H|+1) = 4*|H|^2+|H| group multiplications:
   # loop
   for h1 in H do
      # 2 * |H|  group multiplications
      Hs := h1^H;

      # FOR theta:
      # some elements g has no solution of h1^2=g
      # so We force them to thetas_rec with zero value:
      if not IsBound(thetas_rec.(String(h1))) then
          thetas_rec.(String(h1)) := 0;
      fi;

      # for others collect square solutions counts:
      square_h1:=h1^2;  # 1 group multiplication
      if not IsBound(thetas_rec.(String(square_h1))) then
          thetas_rec.(String(square_h1)) := 1;
      else
          thetas_rec.(String(square_h1)) := thetas_rec.(String(square_h1))+1;
      fi;
      # END FOR
   
      if not IsBound(cents_tmp.(String(h1))) then
          # 2*|H| group multiplications:
          cents_tmp.(String(h1)):=cent(H, h1);
          Print("cents_tmp of ",String(h1), " = ", cents_tmp.(String(h1)),"\n");
      fi;

      # 2th loop
      for h2 in Hs do 
        # Print("h2:", h2,"\n");
        Add(steps, h2); 
      od;
   od;
   Print("steps:", steps,"\n");

   # 5+4+2+4=16/ per step digital multiplications
   # |H|^2*16 in common
   left:=List(steps, h->(3*thetas_rec.(String(h))^4*cents_tmp.(String(h))
                         + 3*thetas_rec.(String(h))^2*cents_tmp.(String(h))^2
                         + cents_tmp.(String(h))^3-3*cents_tmp.(String(h))^4)); # loop
   Print("left: ", left,"\n");
   Print("sum left: ", Sum(left), "\n");
   
   # |H|+2 digital multiplications:
   right:=4*Size(H)*Sum(H, h->cents_tmp.(String(h))^2); # loop

   Print("thetas: ", thetas_rec,"\n");
   Print("cents_rec: ", cents_tmp,"\n");
   Print("right: ", right);
   return(Sum(left)=right);
end;


kolesnikov_SR_sufficient:=function(H)
    return(kolesnikov_left(H)=kolesnikov_right(H));
end;


kolesnikov_th2 := function(H)
        # gap> kolesnikov_th2(ElementaryAbelianGroup(2^2));
        # left: 4096
        # right: 4096
        
        # gap> kolesnikov_th2(Group((1,2)));
        # left: 160
        # right: 160
        
        # gap> kolesnikov_th2(Group((1,2),(1,2,3)));
        # left: 6876
        # right: 6876

        # gap> kolesnikov_th2(Group((1,2),(1,2,3,4)));
        # left: 930240
        # right: 1275840

        # gap> kolesnikov_th2(Group((1,2),(1,2,3,4,5)));
        # left: 196632000
        # right: 645220800

        local idxs, left, right;
        idxs := Flat(List(List(H, elm->elm^H),elm->List(elm)));
        
        left := Sum(List(idxs, elm->
                (3*theta(H, elm)^4*Size(Centralizer(H, elm))
                 + 3*theta(H, elm)^2*Size(Centralizer(H, elm))^2
                 + Size(Centralizer(H, elm))^3)));
        Print("left: ", left, "\n");

        right := (4*Size(H)*Sum(List(H, elm->Size(Centralizer(H, elm))^2))
         + 3 * Sum(List(idxs, elm->Size(Centralizer(H, elm))^4)));
        Print("right: ", right, "\n");
end;

#####################

tests:=function(H, sigma)
    # gap> tests(Group((1,2),(1,2,3)), (1, 4)(2, 5)(3, 6)); 

    local result;
    result:=test_l1(H, sigma);
    Print("test_l1 ", result, "\n");
    result:=test_l2(H, sigma);
    Print("test_l2 ", result, "\n");
    result:=test_l3(H, sigma);
    Print("test_l3 ", result, "\n");
    result:=test_l4(H, sigma);
    Print("\ntest_l4 ", result, "\n");

    result:=(test_sum_5_right(H))=test_sum_5_left(H);
    Print("\ntest_sum_5 ", result, "\n");
    result:=(test_sum_8_right(H))=test_sum_8_left(H);
    Print("\ntest_sum_8 ", result, "\n");
    result:=test_sum_8_1(H);
    Print("\ntest_sum_8_1 ", result, "\n");
    

end;
