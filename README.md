# FastGoat
What C# can do for studying Finite Groups, abelians or not, quotient groups, direct products and many more...

```
var z = new Zn(4, 5);
var g = z.GroupElement(z.CE(1, 0), z.CE(0, 1)).Generate();
var h = z.Monogenic(z.CE(0, 1));
g.SortBy = h.SortBy = SortBy.Value;
g.DisplayElements("G");
h.DisplayElements( "H");
var gh = g.Over(h);
gh.SortBy = SortBy.Value;
gh.DisplayElements();
gh.DisplayGroupTable();
gh.DisplayClasses();
```
Will output
```
|G| = 20 
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0)
a = ( 0, 1)
b = ( 0, 2)
c = ( 0, 3)
d = ( 0, 4)
e = ( 1, 0)
f = ( 1, 1)
g = ( 1, 2)
h = ( 1, 3)
i = ( 1, 4)
j = ( 2, 0)
k = ( 2, 1)
l = ( 2, 2)
m = ( 2, 3)
n = ( 2, 4)
o = ( 3, 0)
p = ( 3, 1)
q = ( 3, 2)
r = ( 3, 3)
s = ( 3, 4)

|H| = 5 
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0)
a = ( 0, 1)
b = ( 0, 2)
c = ( 0, 3)
d = ( 0, 4)

|G/H| = 4 
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0)
a = ( 1, 0)
b = ( 2, 0)
c = ( 3, 0)

|G/H| = 4 
IsGroup      :  True
IsCommutative:  True

 *|@ a b c
--|--------
 @|@ a b c
 a|a b c @
 b|b c @ a
 c|c @ a b

Class of : ( 0, 0)
        ( 0, 0)
        ( 0, 1)
        ( 0, 2)
        ( 0, 3)
        ( 0, 4)
Class of : ( 1, 0)
        ( 1, 0)
        ( 1, 1)
        ( 1, 2)
        ( 1, 3)
        ( 1, 4)
Class of : ( 2, 0)
        ( 2, 0)
        ( 2, 1)
        ( 2, 2)
        ( 2, 3)
        ( 2, 4)
Class of : ( 3, 0)
        ( 3, 0)
        ( 3, 1)
        ( 3, 2)
        ( 3, 3)
        ( 3, 4)
```

Example with S4, A4 and K4
```
var S4 = new Sn(4);
var A4 = S4.GroupElement(S4.C(1, 2, 3), S4.C(2, 3, 4)).Generate();
var Klein = S4.GroupElement(S4.C((1, 2), (3, 4)), S4.C((1, 3), (2, 4))).Generate();
A4.DisplayElements("A4", "in S4");
Klein.DisplayElements("Klein", "in S4");

var Q = A4.Over(Klein);
Q.Details(infos: "in S4");
Q.DisplayClasses();
```

Will output

```
|A4| = 12 in S4
IsGroup      :  True
IsCommutative: False

@[1]  = [ 1 2 3 4](+) Invariants : [1 2 3 4]
a[2]  = [ 2 1 4 3](+) Cycles : (1 2) (3 4)
b[2]  = [ 3 4 1 2](+) Cycles : (1 3) (2 4)
c[2]  = [ 4 3 2 1](+) Cycles : (1 4) (2 3)
d[3]  = [ 1 3 4 2](+) Invariants : [1] Cycles : (2 3 4)
e[3]  = [ 1 4 2 3](+) Invariants : [1] Cycles : (2 4 3)
f[3]  = [ 2 3 1 4](+) Invariants : [4] Cycles : (1 2 3)
g[3]  = [ 2 4 3 1](+) Invariants : [3] Cycles : (1 2 4)
h[3]  = [ 3 1 2 4](+) Invariants : [4] Cycles : (1 3 2)
i[3]  = [ 3 2 4 1](+) Invariants : [2] Cycles : (1 3 4)
j[3]  = [ 4 1 3 2](+) Invariants : [3] Cycles : (1 4 2)
k[3]  = [ 4 2 1 3](+) Invariants : [2] Cycles : (1 4 3)

|Klein| = 4 in S4
IsGroup      :  True
IsCommutative:  True

@[1]  = [ 1 2 3 4](+) Invariants : [1 2 3 4]
a[2]  = [ 2 1 4 3](+) Cycles : (1 2) (3 4)
b[2]  = [ 3 4 1 2](+) Cycles : (1 3) (2 4)
c[2]  = [ 4 3 2 1](+) Cycles : (1 4) (2 3)

|A4/Klein| = 3 in S4
IsGroup      :  True
IsCommutative:  True

@[1]  = [ 1 2 3 4](+) Invariants : [1 2 3 4]
a[3]  = [ 1 3 4 2](+) Invariants : [1] Cycles : (2 3 4)
b[3]  = [ 1 4 2 3](+) Invariants : [1] Cycles : (2 4 3)

|A4/Klein| = 3 in S4
IsGroup      :  True
IsCommutative:  True

 *|@ a b
--|------
 @|@ a b
 a|a b @
 b|b @ a

Class of : [ 1 2 3 4](+) Invariants : [1 2 3 4]
        [ 1 2 3 4](+) Invariants : [1 2 3 4]
        [ 2 1 4 3](+) Cycles : (1 2) (3 4)
        [ 3 4 1 2](+) Cycles : (1 3) (2 4)
        [ 4 3 2 1](+) Cycles : (1 4) (2 3)
Class of : [ 1 3 4 2](+) Invariants : [1] Cycles : (2 3 4)
        [ 1 3 4 2](+) Invariants : [1] Cycles : (2 3 4)
        [ 2 4 3 1](+) Invariants : [3] Cycles : (1 2 4)
        [ 3 1 2 4](+) Invariants : [4] Cycles : (1 3 2)
        [ 4 2 1 3](+) Invariants : [2] Cycles : (1 4 3)
Class of : [ 1 4 2 3](+) Invariants : [1] Cycles : (2 4 3)
        [ 1 4 2 3](+) Invariants : [1] Cycles : (2 4 3)
        [ 2 3 1 4](+) Invariants : [4] Cycles : (1 2 3)
        [ 3 2 4 1](+) Invariants : [2] Cycles : (1 3 4)
        [ 4 1 3 2](+) Invariants : [3] Cycles : (1 4 2)
    
```

Some computing on Z/2Z x Z/2Z x Z/2Z x Z/3Z

```
var z = new Zn(2, 2, 2, 3);
var z24 = z.GenerateAll();
z24.DisplayElements("G");

// Greatest order element of the group
var c6 = z.Monogenic(z.CE(1, 1, 1, 2));
c6.DisplayElements("C6");

// Quotient group 
var k = z24.Over(c6);
k.Details();

// Greatest order element of the quotient group
var c20 = z.Monogenic(z.CE(0, 0, 1, 0));
c20.DisplayElements("C2");

k.Over(c20).Details();

var c21 = z.Monogenic(z.CE(0, 1, 0, 0));
c21.DisplayElements("C2'");

// Direct product of the factors
c20.DirectProduct(c21).DirectProduct(c6).DisplayElements("C2.C2'.C6");
```

Will output

```
|G| = 24 
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 0, 0, 1, 0)
b[2]  = ( 0, 1, 0, 0)
c[2]  = ( 0, 1, 1, 0)
d[2]  = ( 1, 0, 0, 0)
e[2]  = ( 1, 0, 1, 0)
f[2]  = ( 1, 1, 0, 0)
g[2]  = ( 1, 1, 1, 0)
h[3]  = ( 0, 0, 0, 1)
i[3]  = ( 0, 0, 0, 2)
j[6]  = ( 0, 0, 1, 1)
k[6]  = ( 0, 0, 1, 2)
l[6]  = ( 0, 1, 0, 1)
m[6]  = ( 0, 1, 0, 2)
n[6]  = ( 0, 1, 1, 1)
o[6]  = ( 0, 1, 1, 2)
p[6]  = ( 1, 0, 0, 1)
q[6]  = ( 1, 0, 0, 2)
r[6]  = ( 1, 0, 1, 1)
s[6]  = ( 1, 0, 1, 2)
t[6]  = ( 1, 1, 0, 1)
u[6]  = ( 1, 1, 0, 2)
v[6]  = ( 1, 1, 1, 1)
w[6]  = ( 1, 1, 1, 2)

|C6| = 6 
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 1, 1, 1, 0)
b[3]  = ( 0, 0, 0, 1)
c[3]  = ( 0, 0, 0, 2)
d[6]  = ( 1, 1, 1, 1)
e[6]  = ( 1, 1, 1, 2)

|G/C6| = 4 
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 0, 0, 1, 0)
b[2]  = ( 0, 1, 0, 0)
c[2]  = ( 0, 1, 1, 0)

|G/C6| = 4 
IsGroup      :  True
IsCommutative:  True

 *|@ a b c
--|--------
 @|@ a b c
 a|a @ c b
 b|b c @ a
 c|c b a @

|C2| = 2 
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 0, 0, 1, 0)

|G/C6/C2| = 2 
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 0, 1, 0, 0)

|G/C6/C2| = 2 
IsGroup      :  True
IsCommutative:  True

 *|@ a
--|----
 @|@ a
 a|a @

|C2'| = 2 
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 0, 1, 0, 0)

|C2.C2'.C6| = 24 
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 0, 0, 1, 0)
b[2]  = ( 0, 1, 0, 0)
c[2]  = ( 0, 1, 1, 0)
d[2]  = ( 1, 0, 0, 0)
e[2]  = ( 1, 0, 1, 0)
f[2]  = ( 1, 1, 0, 0)
g[2]  = ( 1, 1, 1, 0)
h[3]  = ( 0, 0, 0, 1)
i[3]  = ( 0, 0, 0, 2)
j[6]  = ( 0, 0, 1, 1)
k[6]  = ( 0, 0, 1, 2)
l[6]  = ( 0, 1, 0, 1)
m[6]  = ( 0, 1, 0, 2)
n[6]  = ( 0, 1, 1, 1)
o[6]  = ( 0, 1, 1, 2)
p[6]  = ( 1, 0, 0, 1)
q[6]  = ( 1, 0, 0, 2)
r[6]  = ( 1, 0, 1, 1)
s[6]  = ( 1, 0, 1, 2)
t[6]  = ( 1, 1, 0, 1)
u[6]  = ( 1, 1, 0, 2)
v[6]  = ( 1, 1, 1, 1)
w[6]  = ( 1, 1, 1, 2)

```

Computing Canonical Decomposition of C20 x C30 and C15 x C10 : May the BRUTEFORCE be with you !!!

```
var z20x30 = new Zn(20, 30).GenerateAll();
GroupExt.InvariantsFactors(z20x30);

var z15x20 = new Zn(15, 20).GenerateAll();
GroupExt.InvariantsFactors(z15x20);
```

Will output
```
Invariants factors of G = C20 x C30
|G| = 600 
IsGroup      :  True
IsCommutative:  True

C60 = ( 1, 1); |<C60>|=60
C60 is SubGroup of G : True
|G/C60| = 10 
IsGroup      :  True
IsCommutative:  True

C10 = ( 0, 1); |<C10>|=10
C10 is SubGroup of G/C60 : True
|G/C60/C10| = 1 
IsGroup      :  True
IsCommutative:  True

C20 x C30 = G[600] ~ C10 x C60
-----------------------------
Invariants factors of G = C15 x C20
|G| = 300 
IsGroup      :  True
IsCommutative:  True

C60 = ( 1, 1); |<C60>|=60
C60 is SubGroup of G : True
|G/C60| = 5 
IsGroup      :  True
IsCommutative:  True

C5 = ( 0, 1); |<C5>|=5
C5 is SubGroup of G/C60 : True
|G/C60/C5| = 1 
IsGroup      :  True
IsCommutative:  True

C15 x C20 = G[300] ~ C5 x C60
-----------------------------
```