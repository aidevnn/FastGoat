# FastGoat
What C# can do for studying Finite Groups, abelians or not, quotient groups, direct products and many more...

## New version
This version is faster and it allows to study more complex cases.

## Example
Searching semi-direct product of C3 by C7 in S7

```csharp
var s7 = new Sn(7);
Perm c7 = (s7, (1, 2, 3, 4, 5, 6, 7));
Perm c2 = (s7, (1, 2));
var S7 = Group.Generate(c2, c7);

var allC7 = S7.Where(e => S7.GetOrderOf(e) == 7);
var allC3 = S7.Where(e => S7.GetOrderOf(e) == 3);
var set = allC7.SelectMany(h => allC3.Select(k => (h, k)));
Console.WriteLine("|S7|={0}, |{{HK in S7 with H~C7 and K~C3}}| = {1}", S7.Count(), set.Count());
Console.WriteLine();

var s = set.First(e => Group.Generate(e.h, e.k).Count() == 21);
Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", s.h, s.k);

Console.WriteLine();

var H = Group.Generate(s.h);
var K = Group.Generate(s.k);
var G = Group.Generate(s.h, s.k);
var GoH = G.Over(H);

H.DisplayDetails("H");
K.DisplayDetails("K");
G.DisplayDetails("G=C3⋊C7");
GoH.DisplayDetails("G/H");
GoH.DisplayCosets();
```

will output

```dotnetcli
|S7|=5040, |{HK in S7 with H~C7 and K~C3}| = 252000

First Solution |HK| = 21 : h = [(1 2 3 4 5 6 7)] and k = [(2 3 5)(4 7 6)]

|H| = 7 in S7
is AbelianGroup

Elements
(1)[1] = []
(2)[7] = [(1 2 3 4 5 6 7)]
(3)[7] = [(1 3 5 7 2 4 6)]
(4)[7] = [(1 4 7 3 6 2 5)]
(5)[7] = [(1 5 2 6 3 7 4)]
(6)[7] = [(1 6 4 2 7 5 3)]
(7)[7] = [(1 7 6 5 4 3 2)]

Table
1 2 3 4 5 6 7
2 3 4 5 6 7 1
3 4 5 6 7 1 2
4 5 6 7 1 2 3
5 6 7 1 2 3 4
6 7 1 2 3 4 5
7 1 2 3 4 5 6

|K| = 3 in S7
is AbelianGroup

Elements
(1)[1] = []
(2)[3] = [(2 3 5)(4 7 6)]
(3)[3] = [(2 5 3)(4 6 7)]

Table
1 2 3
2 3 1
3 1 2

|G=C3⋊C7| = 21 in S7
is NotAbelianGroup

Elements
( 1)[ 1] = []
( 2)[ 3] = [(2 3 5)(4 7 6)]
( 3)[ 3] = [(2 5 3)(4 6 7)]
( 4)[ 3] = [(1 2 4)(3 6 5)]
( 5)[ 3] = [(1 2 6)(4 7 5)]
( 6)[ 3] = [(1 3 7)(2 5 4)]
( 7)[ 3] = [(1 3 4)(2 7 6)]
( 8)[ 3] = [(1 4 2)(3 5 6)]
( 9)[ 3] = [(1 4 3)(2 6 7)]
(10)[ 3] = [(1 5 7)(3 6 4)]
(11)[ 3] = [(1 5 6)(2 7 3)]
(12)[ 3] = [(1 6 2)(4 5 7)]
(13)[ 3] = [(1 6 5)(2 3 7)]
(14)[ 3] = [(1 7 5)(3 4 6)]
(15)[ 3] = [(1 7 3)(2 4 5)]
(16)[ 7] = [(1 2 3 4 5 6 7)]
(17)[ 7] = [(1 3 5 7 2 4 6)]
(18)[ 7] = [(1 4 7 3 6 2 5)]
(19)[ 7] = [(1 5 2 6 3 7 4)]
(20)[ 7] = [(1 6 4 2 7 5 3)]
(21)[ 7] = [(1 7 6 5 4 3 2)]

Table
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
 2  3  1  5 16  7 17 18  8 19 10 13 20 15 21  4  6  9 11 12 14
 3  1  2 16  4 17  6  9 18 11 19 20 12 21 14  5  7  8 10 13 15
 4  7 19  8 20 10 21  1 13 16 15  3 17  5 18  6  9 11 12 14  2
 5 17 11 18 12 19 14  2 20  4 21  1  6 16  9  7  8 10 13 15  3
 6 10 16 13 17 15 18 19  3 20  5  7 21  8  1  9 11 12 14  2  4
 7 19  4 20  6 21  9 11  1 12 16 17 14 18  2  8 10 13 15  3  5
 8 21 12  1 14 16  2  4 17  6 18 19  9 20 11 10 13 15  3  5  7
 9 15 20  3 21  5  1 16  7 17  8 10 18 13 19 11 12 14  2  4  6
10 16  6 17  9 18 11 12 19 14 20 21  2  1  4 13 15  3  5  7  8
11  5 17  7 18  8 19 20 10 21 13 15  1  3 16 12 14  2  4  6  9
12  8 21 10  1 13 16 17 15 18  3  5 19  7 20 14  2  4  6  9 11
13 18 14 19  2 20  4  6 21  9  1 16 11 17 12 15  3  5  7  8 10
14 13 18 15 19  3 20 21  5  1  7  8 16 10 17  2  4  6  9 11 12
15 20  9 21 11  1 12 14 16  2 17 18  4 19  6  3  5  7  8 10 13
16  6 10  9 13 11 15  3 12  5 14  2  7  4  8 17 18 19 20 21  1
17 11  5 12  7 14  8 10  2 13  4  6 15  9  3 18 19 20 21  1 16
18 14 13  2 15  4  3  5  6  7  9 11  8 12 10 19 20 21  1 16 17
19  4  7  6  8  9 10 13 11 15 12 14  3  2  5 20 21  1 16 17 18
20  9 15 11  3 12  5  7 14  8  2  4 10  6 13 21  1 16 17 18 19
21 12  8 14 10  2 13 15  4  3  6  9  5 11  7  1 16 17 18 19 20

|G/H| = 3 in S7
is AbelianGroup

Elements
(1)[1] = []
(2)[3] = [(2 3 5)(4 7 6)]
(3)[3] = [(2 5 3)(4 6 7)]

Table
1 2 3
2 3 1
3 1 2

Cosets
(1)[1] = []
      []
      [(1 2 3 4 5 6 7)]
      [(1 3 5 7 2 4 6)]
      [(1 4 7 3 6 2 5)]
      [(1 5 2 6 3 7 4)]
      [(1 6 4 2 7 5 3)]
      [(1 7 6 5 4 3 2)]
(2)[3] = [(2 3 5)(4 7 6)]
      [(2 3 5)(4 7 6)]
      [(1 2 4)(3 6 5)]
      [(1 3 7)(2 5 4)]
      [(1 4 3)(2 6 7)]
      [(1 5 6)(2 7 3)]
      [(1 6 2)(4 5 7)]
      [(1 7 5)(3 4 6)]
(3)[3] = [(2 5 3)(4 6 7)]
      [(2 5 3)(4 6 7)]
      [(1 2 6)(4 7 5)]
      [(1 3 4)(2 7 6)]
      [(1 4 2)(3 5 6)]
      [(1 5 7)(3 6 4)]
      [(1 6 5)(2 3 7)]
      [(1 7 3)(2 4 5)]

```

## References

<b>Daniel Guin, Thomas Hausberger.</b> ALGÈBRE T1 Groupes, corps et théorie de Galois. 