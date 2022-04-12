# FastGoat
Fastest C# can do for studying Finite Groups, abelians or not, quotient groups, direct products and many more...

```
var z = new Zn(4, 5);
var g = z.Union(z.CreateElement(1, 0),z.CreateElement(0, 1)).Develop("G");
var h = z.Monogenic(z.CreateElement(0, 1), "H");
var gh = g.Over(h);
g.SortBy = h.SortBy = gh.SortBy = SortBy.Value;
g.DisplayElements();
h.DisplayElements();
gh.Details();
gh.DisplayClasses();
```
Will output
```
|G| = 20 in Z/4Z x Z/5Z
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0)[1]
a = ( 0, 1)[5]
b = ( 0, 2)[5]
c = ( 0, 3)[5]
d = ( 0, 4)[5]
e = ( 1, 0)[4]
f = ( 1, 1)[20]
g = ( 1, 2)[20]
h = ( 1, 3)[20]
i = ( 1, 4)[20]
j = ( 2, 0)[2]
k = ( 2, 1)[10]
l = ( 2, 2)[10]
m = ( 2, 3)[10]
n = ( 2, 4)[10]
o = ( 3, 0)[4]
p = ( 3, 1)[20]
q = ( 3, 2)[20]
r = ( 3, 3)[20]
s = ( 3, 4)[20]

|H| = 5 in Z/4Z x Z/5Z
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0)[1]
a = ( 0, 1)[5]
b = ( 0, 2)[5]
c = ( 0, 3)[5]
d = ( 0, 4)[5]

|G/H| = 4 with |G| = 20 and |H| = 5, OpBoth
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0)[1]
a = ( 1, 0)[4]
b = ( 2, 0)[2]
c = ( 3, 0)[4]

|G/H| = 4 with |G| = 20 and |H| = 5, OpBoth
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
var A4 = S4.Union(S4.Cycle(1, 2, 3), S4.Cycle(2, 3, 4)).Develop("A4");
var Klein = S4.Union(S4.Cycles((1, 2), (3, 4)), S4.Cycles((1, 3), (2, 4))).Develop("Klein");
var Q = A4.Over(Klein);
A4.DisplayElements();
Klein.DisplayElements();
Q.Details();
Q.DisplayClasses();
```

Will output

```
|A4| = 12 in S4
IsGroup      :  True
IsCommutative: False

@ = ( 1  2  3  4)[1](+) Invariants : [1 2 3 4]
a = ( 2  1  4  3)[2](+) Cycles : (1 2) (3 4)
b = ( 3  4  1  2)[2](+) Cycles : (1 3) (2 4)
c = ( 4  3  2  1)[2](+) Cycles : (1 4) (2 3)
d = ( 1  3  4  2)[3](+) Invariants : [1] Cycles : (2 3 4)
e = ( 1  4  2  3)[3](+) Invariants : [1] Cycles : (2 4 3)
f = ( 2  3  1  4)[3](+) Invariants : [4] Cycles : (1 2 3)
g = ( 2  4  3  1)[3](+) Invariants : [3] Cycles : (1 2 4)
h = ( 3  1  2  4)[3](+) Invariants : [4] Cycles : (1 3 2)
i = ( 3  2  4  1)[3](+) Invariants : [2] Cycles : (1 3 4)
j = ( 4  1  3  2)[3](+) Invariants : [3] Cycles : (1 4 2)
k = ( 4  2  1  3)[3](+) Invariants : [2] Cycles : (1 4 3)

|Klein| = 4 in S4
IsGroup      :  True
IsCommutative:  True

@ = ( 1  2  3  4)[1](+) Invariants : [1 2 3 4]
a = ( 2  1  4  3)[2](+) Cycles : (1 2) (3 4)
b = ( 3  4  1  2)[2](+) Cycles : (1 3) (2 4)
c = ( 4  3  2  1)[2](+) Cycles : (1 4) (2 3)

|A4/Klein| = 3 with |A4| = 12 and |Klein| = 4, OpBoth
IsGroup      :  True
IsCommutative:  True

@ = ( 1  2  3  4)[1](+) Invariants : [1 2 3 4]
a = ( 1  3  4  2)[3](+) Invariants : [1] Cycles : (2 3 4)
b = ( 1  4  2  3)[3](+) Invariants : [1] Cycles : (2 4 3)

|A4/Klein| = 3 with |A4| = 12 and |Klein| = 4, OpBoth
 *|@ a b
--|------
 @|@ a b
 a|a b @
 b|b @ a


Class of : ( 1  2  3  4)
    ( 1  2  3  4)
    ( 2  1  4  3)
    ( 3  4  1  2)
    ( 4  3  2  1)
Class of : ( 1  3  4  2)
    ( 1  3  4  2)
    ( 2  4  3  1)
    ( 3  1  2  4)
    ( 4  2  1  3)
Class of : ( 1  4  2  3)
    ( 1  4  2  3)
    ( 2  3  1  4)
    ( 3  2  4  1)
    ( 4  1  3  2)
    
```

Invariants factors of Z/2Z x Z/2Z x Z/2Z x Z/3Z

```
var z = new Zn(2, 2, 2, 3);
var z24 = z.Union(z.BaseCanonic).Develop("G");
z24.DisplayElements();

// Greatest order element of the group
var c6 = z.Monogenic(z.CE(1, 1, 1, 2), "C6");
c6.DisplayElements();

// Quotient group 
var k = z24.Over(c6);
k.Details();

// Greatest order element of the quotient group
var c20 = z.Monogenic(z.CE(0, 0, 1, 0), "C2");
k.Over(c20).Details();

var c21 = z.Monogenic(z.CE(0, 1, 0, 0), "C2'");

// Direct product of the invariants factors
c21.DirectProduct(c20).DirectProduct(c6).DisplayElements();
```

Will output

```
|G| = 24 in Z/2Z x Z/2Z x Z/2Z x Z/3Z
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0, 0, 0)[1]
a = ( 0, 0, 1, 0)[2]
b = ( 0, 1, 0, 0)[2]
c = ( 0, 1, 1, 0)[2]
d = ( 1, 0, 0, 0)[2]
e = ( 1, 0, 1, 0)[2]
f = ( 1, 1, 0, 0)[2]
g = ( 1, 1, 1, 0)[2]
h = ( 0, 0, 0, 1)[3]
i = ( 0, 0, 0, 2)[3]
j = ( 0, 0, 1, 1)[6]
k = ( 0, 0, 1, 2)[6]
l = ( 0, 1, 0, 1)[6]
m = ( 0, 1, 0, 2)[6]
n = ( 0, 1, 1, 1)[6]
o = ( 0, 1, 1, 2)[6]
p = ( 1, 0, 0, 1)[6]
q = ( 1, 0, 0, 2)[6]
r = ( 1, 0, 1, 1)[6]
s = ( 1, 0, 1, 2)[6]
t = ( 1, 1, 0, 1)[6]
u = ( 1, 1, 0, 2)[6]
v = ( 1, 1, 1, 1)[6]
w = ( 1, 1, 1, 2)[6]

|C6| = 6 in Z/2Z x Z/2Z x Z/2Z x Z/3Z
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0, 0, 0)[1]
a = ( 1, 1, 1, 0)[2]
b = ( 0, 0, 0, 1)[3]
c = ( 0, 0, 0, 2)[3]
d = ( 1, 1, 1, 1)[6]
e = ( 1, 1, 1, 2)[6]

|G/C6| = 4 with |G| = 24 and |C6| = 6
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0, 0, 0)[1]
a = ( 0, 0, 1, 0)[2]
b = ( 0, 1, 0, 0)[2]
c = ( 0, 1, 1, 0)[2]

|G/C6| = 4 with |G| = 24 and |C6| = 6
 *|@ a b c
--|--------
 @|@ a b c
 a|a @ c b
 b|b c @ a
 c|c b a @


|G/C6/C2| = 2 with |G/C6| = 4 and |C2| = 2
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0, 0, 0)[1]
a = ( 0, 1, 0, 0)[2]

|G/C6/C2| = 2 with |G/C6| = 4 and |C2| = 2
 *|@ a
--|----
 @|@ a
 a|a @


|C2'.C2.C6| = 24 in Z/2Z x Z/2Z x Z/2Z x Z/3Z
IsGroup      :  True
IsCommutative:  True

@ = ( 0, 0, 0, 0)[1]
a = ( 0, 0, 1, 0)[2]
b = ( 0, 1, 0, 0)[2]
c = ( 0, 1, 1, 0)[2]
d = ( 1, 0, 0, 0)[2]
e = ( 1, 0, 1, 0)[2]
f = ( 1, 1, 0, 0)[2]
g = ( 1, 1, 1, 0)[2]
h = ( 0, 0, 0, 1)[3]
i = ( 0, 0, 0, 2)[3]
j = ( 0, 0, 1, 1)[6]
k = ( 0, 0, 1, 2)[6]
l = ( 0, 1, 0, 1)[6]
m = ( 0, 1, 0, 2)[6]
n = ( 0, 1, 1, 1)[6]
o = ( 0, 1, 1, 2)[6]
p = ( 1, 0, 0, 1)[6]
q = ( 1, 0, 0, 2)[6]
r = ( 1, 0, 1, 1)[6]
s = ( 1, 0, 1, 2)[6]
t = ( 1, 1, 0, 1)[6]
u = ( 1, 1, 0, 2)[6]
v = ( 1, 1, 1, 1)[6]


```