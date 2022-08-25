# FastGoat
What C# can do for studying Finite Groups, abelians or not, quotient groups, direct products and many more...

### Starting with Abelian Group
```
var z = new Zn(4, 5);
var g = z.GroupElement(z.CE(1, 0), z.CE(0, 1)).Generate();
var h = g.Monogenic(z.CE(0, 1));
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

### Example with S4, A4 and K4
```
var S4 = new Sn(4);
var A4 = S4.GroupElement(S4.C(1, 2, 3), S4.C(2, 3, 4)).Generate();
var Klein = A4.GroupElement(S4.C((1, 2), (3, 4)), S4.C((1, 3), (2, 4))).Generate();
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

### Some computing on Z/2Z x Z/2Z x Z/2Z x Z/3Z

```
var z = new Zn(2, 2, 2, 3);
var z24 = z.GenerateAll();
z24.DisplayElements("G", "Cartesian product Z/2Z x Z/2Z x Z/2Z x Z/3Z");

// Greatest order element of the group
var c6 = z24.Monogenic(z.CE(1, 1, 1, 2));
c6.DisplayElements("C6");

// Quotient group 
var q0 = z24.Over(c6);
q0.Details();

// Greatest order element of the quotient group
var c20 = q0.Monogenic(z.CE(0, 0, 1, 0));
c20.DisplayElements("C2");

var q1 = q0.Over(c20);
q1.Details();

var c21 = q1.Monogenic(z.CE(0, 1, 0, 0));
c21.DisplayElements("C2'");

// Direct product of the invariants factors
Console.WriteLine("###########");
c6.DirectProduct(c20).DisplayElements("C6.C2");
Console.WriteLine("###########");
c6.DirectProduct(c20).DirectProduct(c21).DisplayElements("C6.C2.C2'");
```

Will output

```
|G| = 24 Cartesian product Z/2Z x Z/2Z x Z/2Z x Z/3Z
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

###########
|C6.C2| = 12
IsGroup      :  True
IsCommutative:  True

@[1]  = ( 0, 0, 0, 0)
a[2]  = ( 0, 0, 1, 0)
b[2]  = ( 1, 1, 0, 0)
c[2]  = ( 1, 1, 1, 0)
d[3]  = ( 0, 0, 0, 1)
e[3]  = ( 0, 0, 0, 2)
f[6]  = ( 0, 0, 1, 1)
g[6]  = ( 0, 0, 1, 2)
h[6]  = ( 1, 1, 0, 1)
i[6]  = ( 1, 1, 0, 2)
j[6]  = ( 1, 1, 1, 1)
k[6]  = ( 1, 1, 1, 2)

###########
|C6.C2.C2'| = 24
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

### Computing Invariant Factor Decomposition

```
var z14x21 = Zn.CartesianProduct(14, 21);
GroupExt.InvariantFactors(z14x21);

var z20x30 = Zn.CartesianProduct(20, 30);
GroupExt.InvariantFactors(z20x30);

var z8x18x30 = Zn.CartesianProduct(8, 18, 30); // May the BRUTEFORCE be with you
GroupExt.InvariantFactors(z8x18x30);

var s6 = new Sn(6);
var H0 = s6.GroupElement(s6.KCycle(2), s6.KCycle(3, 4)).Generate("H0");
H0.DisplayElements();
GroupExt.InvariantFactors(H0);

var s9 = new Sn(9);
var H1 = s9.GroupElement(s9.KCycle(3), s9.KCycle(4, 6)).Generate("H1");
H1.DisplayElements();
GroupExt.InvariantFactors(H1);
```

Will output
```
Invariants factors of G = C14 x C21
|G| = 294 
IsGroup      :  True
IsCommutative:  True

C42 = ( 1, 1); |<C42>|=42
C42 is SubGroup of G : True
|G/C42| = 7 
IsGroup      :  True
IsCommutative:  True

C7 = ( 0, 1); |<C7>|=7
C7 is SubGroup of G/C42 : True
|G/C42/C7| = 1 
IsGroup      :  True
IsCommutative:  True

C14 x C21 = G[294] ~ C7 x C42
-----------------------------
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
Invariants factors of G = C8 x C18 x C30
|G| = 4320 
IsGroup      :  True
IsCommutative:  True

C360 = ( 1, 1, 1); |<C360>|=360
C360 is SubGroup of G : True
|G/C360| = 12 
IsGroup      :  True
IsCommutative:  True

C6 = ( 0, 0, 1); |<C6>|=6
C6 is SubGroup of G/C360 : True
|G/C360/C6| = 2 
IsGroup      :  True
IsCommutative:  True

C2 = ( 0, 1, 0); |<C2>|=2
C2 is SubGroup of G/C360/C6 : True
|G/C360/C6/C2| = 1 
IsGroup      :  True
IsCommutative:  True

C8 x C18 x C30 = G[4320] ~ C2 x C6 x C360
-----------------------------
|H0| = 8 
IsGroup      :  True
IsCommutative:  True

@[1]  = [ 1 2 3 4 5 6](+) Invariants : [1 2 3 4 5 6]
a[2]  = [ 1 2 5 6 3 4](+) Invariants : [1 2] Cycles : (3 5) (4 6)
b[2]  = [ 2 1 3 4 5 6](-) Invariants : [3 4 5 6] Cycles : (1 2)
c[2]  = [ 2 1 5 6 3 4](-) Cycles : (1 2) (3 5) (4 6)
d[4]  = [ 1 2 4 5 6 3](-) Invariants : [1 2] Cycles : (3 4 5 6)
e[4]  = [ 1 2 6 3 4 5](-) Invariants : [1 2] Cycles : (3 6 5 4)
f[4]  = [ 2 1 4 5 6 3](+) Cycles : (1 2) (3 4 5 6)
g[4]  = [ 2 1 6 3 4 5](+) Cycles : (1 2) (3 6 5 4)

Invariants factors of G = H0
|G| = 8 
IsGroup      :  True
IsCommutative:  True

C4 = [ 1 2 4 5 6 3](-) Invariants : [1 2] Cycles : (3 4 5 6); |<C4>|=4
C4 is SubGroup of G : True
|G/C4| = 2 
IsGroup      :  True
IsCommutative:  True

C2 = [ 2 1 3 4 5 6](-) Invariants : [3 4 5 6] Cycles : (1 2); |<C2>|=2
C2 is SubGroup of G/C4 : True
|G/C4/C2| = 1 
IsGroup      :  True
IsCommutative:  True

H0 = G[8] ~ C2 x C4
-----------------------------
|H1| = 18 
IsGroup      :  True
IsCommutative:  True

@[1]  = [ 1 2 3 4 5 6 7 8 9](+) Invariants : [1 2 3 4 5 6 7 8 9]
a[2]  = [ 1 2 3 7 8 9 4 5 6](-) Invariants : [1 2 3] Cycles : (4 7) (5 8) (6 9)
b[3]  = [ 1 2 3 6 7 8 9 4 5](+) Invariants : [1 2 3] Cycles : (4 6 8) (5 7 9)
c[3]  = [ 1 2 3 8 9 4 5 6 7](+) Invariants : [1 2 3] Cycles : (4 8 6) (5 9 7)
d[3]  = [ 2 3 1 4 5 6 7 8 9](+) Invariants : [4 5 6 7 8 9] Cycles : (1 2 3)
e[3]  = [ 2 3 1 6 7 8 9 4 5](+) Cycles : (1 2 3) (4 6 8) (5 7 9)
f[3]  = [ 2 3 1 8 9 4 5 6 7](+) Cycles : (1 2 3) (4 8 6) (5 9 7)
g[3]  = [ 3 1 2 4 5 6 7 8 9](+) Invariants : [4 5 6 7 8 9] Cycles : (1 3 2)
h[3]  = [ 3 1 2 6 7 8 9 4 5](+) Cycles : (1 3 2) (4 6 8) (5 7 9)
i[3]  = [ 3 1 2 8 9 4 5 6 7](+) Cycles : (1 3 2) (4 8 6) (5 9 7)
j[6]  = [ 1 2 3 5 6 7 8 9 4](-) Invariants : [1 2 3] Cycles : (4 5 6 7 8 9)
k[6]  = [ 1 2 3 9 4 5 6 7 8](-) Invariants : [1 2 3] Cycles : (4 9 8 7 6 5)
l[6]  = [ 2 3 1 5 6 7 8 9 4](-) Cycles : (1 2 3) (4 5 6 7 8 9)
m[6]  = [ 2 3 1 7 8 9 4 5 6](-) Cycles : (1 2 3) (4 7) (5 8) (6 9)
n[6]  = [ 2 3 1 9 4 5 6 7 8](-) Cycles : (1 2 3) (4 9 8 7 6 5)
o[6]  = [ 3 1 2 5 6 7 8 9 4](-) Cycles : (1 3 2) (4 5 6 7 8 9)
p[6]  = [ 3 1 2 7 8 9 4 5 6](-) Cycles : (1 3 2) (4 7) (5 8) (6 9)
q[6]  = [ 3 1 2 9 4 5 6 7 8](-) Cycles : (1 3 2) (4 9 8 7 6 5)

Invariants factors of G = H1
|G| = 18 
IsGroup      :  True
IsCommutative:  True

C6 = [ 1 2 3 5 6 7 8 9 4](-) Invariants : [1 2 3] Cycles : (4 5 6 7 8 9); |<C6>|=6
C6 is SubGroup of G : True
|G/C6| = 3 
IsGroup      :  True
IsCommutative:  True

C3 = [ 2 3 1 4 5 6 7 8 9](+) Invariants : [4 5 6 7 8 9] Cycles : (1 2 3); |<C3>|=3
C3 is SubGroup of G/C6 : True
|G/C6/C3| = 1 
IsGroup      :  True
IsCommutative:  True

H1 = G[18] ~ C3 x C6
-----------------------------

```