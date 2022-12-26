# FastGoat
What C# can do for studying Finite Groups, abelians or not, quotient groups, cartesian products, direct products, semi-direct products and many more...

## Example
Searching semi-direct product $C_7 \rtimes C_3$ in $\textbf{S}_7$

```csharp
GlobalStopWatch.Restart();
var s7 = new Sn(7);
var a = s7[(1, 2, 3, 4, 5, 6, 7)];
var allC3 = s7.Where(p => (p ^ 3) == s7.Neutral()).ToArray();
var b = allC3.First(p => Group.GenerateElements(s7, a, p).Count() == 21);
GlobalStopWatch.Stop();

Console.WriteLine("|S7|={0}, |{{b in S7 with b^3 = 1}}| = {1}",s7.Count(), allC3.Count());
Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", a, b);
Console.WriteLine();

var h = Group.Generate("H", a);
var g21 = Group.Generate("G21", a, b);
DisplayGroup.Head(g21);
DisplayGroup.Head(g21.Over(h));
GlobalStopWatch.Show("Group21");
```

will output

```dotnetcli
|S7|=5040, |{b in S7 with b^3 = 1}| = 351
First Solution |HK| = 21 : h = [(1 2 3 4 5 6 7)] and k = [(2 3 5)(4 7 6)]

|G21| = 21
Type        NonAbelianGroup
BaseGroup   S7

|G21/H| = 3
Type        AbelianGroup
BaseGroup   S7
SuperGroup  |G21| = 21
NormalGroup |H| = 7

# Group21 Time:65 ms
```

## Another Example
Comparing the previous results with the group presented by $\langle (a,\ b) \ | \ a^7=b^3=1,\ a^2=bab^{-1} \rangle$

```csharp
GlobalStopWatch.Restart();
var wg = new WordGroup("a7, b3, a2 = bab-1");
GlobalStopWatch.Stop();

DisplayGroup.Head(wg);
var n = Group.Generate("<a>", wg, wg["a"]);
DisplayGroup.Head(wg.Over(n));
GlobalStopWatch.Show($"{wg}");
Console.WriteLine();
```

will produce

```dotnetcli
|WG[a,b]| = 21
Type        NonAbelianGroup
BaseGroup   WG[a,b]

|WG[a,b]/<a>| = 3
Type        AbelianGroup
BaseGroup   WG[a,b]/<a>
Group           |WG[a,b]| = 21
NormalSubGroup  |<a>| = 7

# WG[a,b] Time:20 ms
```

## Semidirect product using group action

Another way for the previous example
```csharp
GlobalStopWatch.Restart();
var c7 = new Cn(7);
var c3 = new Cn(3);
var g21 = Group.SemiDirectProd(c7, c3);
GlobalStopWatch.Stop();

var n = Group.Generate("N", g21, g21[1, 0]);
DisplayGroup.HeadSdp(g21);
DisplayGroup.Head(g21.Over(n));
GlobalStopWatch.Show("Group21");
```
will output
```dotnetcli
|C7 x: C3| = 21
Type        NonAbelianGroup
BaseGroup    Z7 x Z3
NormalGroup  |C7| = 7
ActionGroup  |C3| = 3

Actions
g=0 y(g) = (0->0, 1->1, 2->2, 3->3, 4->4, 5->5, 6->6)
g=1 y(g) = (0->0, 1->2, 2->4, 3->6, 4->1, 5->3, 6->5)
g=2 y(g) = (0->0, 1->4, 2->1, 3->5, 4->2, 5->6, 6->3)

|(C7 x: C3)/N| = 3
Type        AbelianGroup
BaseGroup   Z7 x Z3
SuperGroup  |C7 x: C3| = 21
NormalGroup |N| = 7

# Group21 Time:0 ms
```

## References

[ALGÈBRE T1](https://laboutique.edpsciences.fr/produit/63/9782759803316/)
<b>Daniel Guin, Thomas Hausberger.</b>
ALGÈBRE T1 Groupes, corps et théorie de Galois.
EDP Sciences.
<i>Premiere Partie.</i>

[Algebra (3rd ed.)](https://bookstore.ams.org/chel-330/)
<b>Saunders MacLane, Garrett Birkhoff.</b>
Algebra (3rd ed.).
American Mathematical Society.
<i>Chapter XII.2 Groups extensions.</i>

[A Course on Finite Groups](https://link.springer.com/book/10.1007/978-1-84882-889-6)
<b>H.E. Rose.</b>
A Course on Finite Groups.
Springer.
<i>All Chapters.</i>

[GroupNames](https://people.maths.bris.ac.uk/~matyd/GroupNames/index.html)
<b>Tim Dokchitser.</b>
Group Names.
<i>Beta</i>

[GAP2022](https://www.gap-system.org)
<b>The GAP Group.</b>
GAP -- Groups, Algorithms, and Programming.
<i>Version 4.12.0; 2022.</i>

[Paper (pdf)](http://www.math.cornell.edu/~kbrown/7350/toddcox.pdf)
<b>Ken Brown.</b>
Mathematics 7350.
Cornell University.
<i>The Todd–Coxeter procedure</i>

[Conway polynomials](http://www.math.rwth-aachen.de/~Frank.Luebeck/data/ConwayPol/index.html)
<b>Frank Lübeck.</b>
Conway polynomials for finite fields
<i>Online data</i>

[AECF](https://hal.archives-ouvertes.fr/AECF/)
<b>Alin Bostan, Frédéric Chyzak, Marc Giusti, Romain Lebreton, Grégoire Lecerf, Bruno Salvy, Éric Schost.</b>
Algorithmes Efficaces en Calcul Formel.
Édition web 1.1
<i>Chapitre IV Factorisation des polynômes.</i>

[hal-01444183](https://hal.archives-ouvertes.fr/hal-01444183)
<b>Xavier Caruso</b>
Computations with p-adic numbers. Journées Nationales de Calcul Formel, In press,
Les cours du CIRM.
<i>Several implementations of p-adic numbers</i>

[Algebraic Factoring](https://dl.acm.org/doi/10.1145/800205.806338)
<b>Barry Trager</b>
Algebraic Factoring and Rational Function Integration
Laboratory for Computer Science, HIT
<i>Norms and Algebraic Factoring</i>