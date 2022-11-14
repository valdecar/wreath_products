### wreath_products

This project contains some functions to check that the finite group $G$ is simply reducible i.e. that two conditions are hold [1]:

1. $\forall g \in G \ \exist \sigma \in G: \ g^{\sigma}=g^{-1}$

2. Kronecker product of every irreducible representation contains each irreducible representation no more than once. 

in other words second condition means that in equation:

$D^{\mu}\times D^{\nu} = \sum_{\sigma}c_{\mu \nu \sigma}D^{\sigma}$

all $c_{\mu \nu \sigma}$ either 0 or 1.

##### E.P. Wigner criteria
According to E.P. Wigner, 1943 [2] this conditions means that for a group $G$ hold:

$\sum_{g \in G}(\theta_{G}^{3}(g)-|C_{G}(g)|^{2})=0$

where $\theta_{G}(g)=|\{h\in G| h^{2}=g\}|$.

This condition is checked by the GAP function `is_SR2` of a file `sr.g`, to use
it run:
```
gap> Read("sr.g");
gap> g:=Group((1,2),(1,2,3));;
gap> is_SR2(g);;
```
which has
```
# complexity:
# |G|*(2*|G| + 1) = 2*|G|^2+|G| group multiplications
# |G|+2*|G|= 3*|G| digital multiplications
```

##### Character criteria
Condition 2 of definition of SR group can be checket using orthogonality of groups irreducible characters (see (5.42) in [1]):

$a_{k}=\frac{1}{|G|}\sum_{s}g_{s} \tilde{\chi_{k}}(\chi_{i}\times\chi_{j})$

This condition is checked by the GAP function `is_SR` of a file `is_SR`, to use it run:
```
gap> Read("sr.g");
gap> g:=Group((1,2),(1,2,3));;
gap> is_SR(g);;
```

##### Wreath product

According to the first theorem of Kolesnikov S.G [3] the wreath product $G=H\wr K$ can be simply reducible only if group $K$ is elementary Abelian 2-group i.e. $K \cong S_{2}$. So, in other words, it is necessary for wreath product to be SR group only if it has the form $H\wr S_{2}$. So it entirely depends on component group $H$ if $H\wr S_{2}$ belongs to SR group or not.

According to the second theorem of Kolesnikov S.G [3] $H\wr S_{2}$ is an SR-group if and only if group $H$ is SR group and it hold a property:

$\sum_{(u,v)\in H\times H, u\in v^{H}} (3\theta_{H}^{4}(u)|C_{H}(u)|+3\theta_{H}^2(u)|C_{H}(u)|^2+|C_{H}(u)|^3)$

=

$4 \|H\|  \sum_{h \in H}|C_{H}(h)|^{2}+3\sum_{(h,f)\in H\times H, h \in f^{H}}|C_{H}(h)|^{4}$

This property is checked by the `GAP` function `kolesnikov_SR_sufficient1` in a file `tests_th2.g` to use it run:
```
gap> Read("tests_th2.g");
gap> d4:=DihedralGroup(4);;
gap> is_SR(d4);; # or is_SR2(d4)
gap> kolesnikov_SR_sufficient1(d4);;
```
Do not forget that this only prove that $D_{4}\wr S_{2}$ is SR group because $D_{4}$ is also such, and because `kolesnikov_SR_sufficient1` condition is hold.

```
# complexity:
# 4*|H|^2+|H| group multiplications
# 16*|H|^2 + |H| + 2 digital multiplications
```

##### Others
There is also some useful functions in a file `sketch.txt` one can use for representations of groups.
There is some useful methods to build the representation of a group in the `SageMath`  in a file `notebooks/representation.ipynb`. 

##### References:

1. M. Hamermesh, Group Theory and its Application to Physical Problems, Dover Books on Physics and Chemistry, Reprint edition.

2. E.P Wigner, On Representations of Certain Finite Groups, American Journal of Mathematics 63, 57-63 (1943).

3. S.G. Kolesnikov, On necessary and sufficient conditions of simply reducibility of wreath product of finite groups, Siberian State University of Science and Thechnology, vol. 19, number 2, 2018, pp. 212-216.