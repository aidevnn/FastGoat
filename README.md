# FastGoat
What C# can do for studying Finite Groups, abelians or not, quotient groups, cartesian products, direct products, semi-direct products and many more...

## Example
Searching semi-direct product $C_7 \rtimes C_3$ in $\textbf{S}_7$

```csharp
GlobalStopWatch.Restart();
var s7 = new Sn(7); // Trivial s7
var a = s7[(1, 2, 3, 4, 5, 6, 7)];
var S7 = Group.GenerateElements(s7, s7[(1, 2)], a);
var allC3 = S7.Where(p => (p ^ 3) == s7.Neutral()).ToArray();
var b = allC3.First(p => Group.GenerateElements(s7, a, p).Count() == 21);
GlobalStopWatch.Stop();

Console.WriteLine("|S7|={0}, |{{b in S7 with b^3 = 1}}| = {1}",S7.Count(), allC3.Count());
Console.WriteLine("First Solution |HK| = 21 : h = {0} and k = {1}", a, b);
Console.WriteLine();

var h = Group.Generate("H", a);
var g21 = Group.Generate(a, b);
DisplayGroup.Head(g21);
DisplayGroup.Head(g21.Over(h));
GlobalStopWatch.Show("Group21");
```

will output

```dotnetcli
|S7|=5040, |{b in S7 with b^3 = 1}| = 351
First Solution |HK| = 21 : h = [(1 2 3 4 5 6 7)] and k = [(1 2 4)(3 6 5)]

|G| = 21
Type        NonAbelianGroup
BaseGroup   S7

|G/H| = 3
Type        AbelianGroup
BaseGroup   S7
SuperGroup  |G| = 21
NormalGroup |H| = 7

# Group21 Time:322 ms
```

## Another Example
Comparing the previous results with the group presented by $\langle (a,\ b) \ | \ a^7=b^3=1,\ a^2=bab^{-1} \rangle$

```csharp

var b = allC3.First(p => (a ^ 2) == p * a * (p ^ -1));
Console.WriteLine("First Solution a^7 = b^3 = 1 and a^2 = b * a * b^-1 : a = {0} and b = {1}", a, b);
```

will output

```dotnetcli
|S7|=5040, |{b in S7 with b^3 = 1}| = 351
First Solution a^7 = b^3 = 1 and a^2 = b * a * b^-1 : a = [(1 2 3 4 5 6 7)] and b = [(1 2 6)(4 7 5)]

|G| = 21
Type        NonAbelianGroup
BaseGroup   S7

|G/H| = 3
Type        AbelianGroup
BaseGroup   S7
SuperGroup  |G| = 21
NormalGroup |H| = 7

# Group21 Time:18 ms
```

## Semidirect product using group action

Another way for the previous example
```csharp
GlobalStopWatch.Restart();
var z7 = new Zn(7); // Trivial Z7
var z3 = new Zn(3); // Trivial Z3
var c7 = Group.Generate("C7", z7[1]);
var c3 = Group.Generate("C3", z3[1]);
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
g=0 y(g):x->x
g=1 y(g):x->x^2
g=2 y(g):x->x^4

|(C7 x: C3)/N| = 3
Type        AbelianGroup
BaseGroup   Z7 x Z3
SuperGroup  |C7 x: C3| = 21
NormalGroup |N| = 7

# Group21 Time:0 ms
```

## References

<b>Daniel Guin, Thomas Hausberger.</b> ALGÈBRE T1 Groupes, corps et théorie de Galois. 
<b>Saunders MacLane, Garrett Birkhoff.<b> Algebra (3rd ed.). American Mathematical Society.