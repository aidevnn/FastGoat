using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection.Emit;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

int[][] AbSubTypes(int[] type)
{
    var all = type.Select(t => Dividors(t).Append(t).ToArray()).MultiLoop()
        .Select(l => l.Order())
        .ToHashSet(new SequenceEquality<int>())
        .Select(l => l.ToArray())
        .Select(l => l.Where(e => e != 1).ToArray())
        .Where(l => l.Length != 0)
        .Append([1])
        .OrderBy(l => l.Length)
        .ThenBy(l => l, Comparer<int[]>.Create((l0, l1) => l0.SequenceCompareTo(l1)))
        .ToArray();

    return all;
}

void example()
{
    var E = new EllGroup<Rational>("-36", "0");
    var O = E.O;
    EllPt<Rational> P = ("-3", "9");
    EllPt<Rational> Q = ("-2", "8");
    Console.WriteLine(new { O, P, Q });
    Console.WriteLine($"-P = {E.Invert(P)}");
    Console.WriteLine($"P + Q = {E.Op(P, Q)}");
    Console.WriteLine($"2P = {E.Times(P, 2)}");
    Console.WriteLine($"2Q = {E.Times(Q, 2)}");
    Console.WriteLine($"2P + 2Q = {E.Op(E.Times(P, 2), E.Times(Q, 2))}");
    Console.WriteLine($"2(P + Q) = {E.Times(E.Op(P, Q), 2)}");
}

int[] EllFp(int a, int b, int p, bool show = false)
{
    var (A, B) = (new ZnInt(p, a), new ZnInt(p, b));
    var disc = 4 * A.Pow(3) + 27 * B.Pow(2);
    if (disc.IsZero())
        return [];

    var E = new EllGroup<ZnInt>(A, B);
    var ell = p.Range().Select(k => new ZnInt(p, k)).ToArray().Grid2D()
        .Select(e => new EllPt<ZnInt>(e.t1, e.t2))
        .Where(e => E.Contains(e.X, e.Y))
        .Distinct()
        .Order()
        .ToArray();

    var gEll = Group.Generate(E, ell);
    if (show)
        DisplayGroup.HeadElements(gEll);

    var abType = Group.AbelianGroupType(gEll);
    Console.WriteLine($"{gEll} ~ {abType.Glue(" x ", "C{0}")}");
    if (gEll.Any(e => !e.IsO && !E.Contains(e.X, e.Y)))
        throw new();

    return abType;
}

void Ellmorph(int a, int b, int nbPrimes = 10, bool show = false)
{
    var disc = 4 * a.Pow(3) + 27 * b.Pow(2);
    var allTypes = Primes10000.Where(p => (2 * disc) % p != 0).Take(nbPrimes)
        .Select(n => EllFp(a, b, n, show))
        .ToHashSet(new SequenceEquality<int>())
        .ToArray();

    var allSubTypes = allTypes.Select(l => AbSubTypes(l.ToArray()).ToHashSet(new SequenceEquality<int>())).ToArray();
    var set = new HashSet<IEnumerable<int>>(new SequenceEquality<int>());
    foreach (var sub in allSubTypes)
    {
        if (set.Count == 0)
        {
            set = sub;
            continue;
        }

        set.IntersectWith(sub);
    }

    allTypes.Println(e => e.Glue(" x ", "C{0}"), $"E[{a},{b}](Rational) ->");
    set.Select(l => l.ToArray()).Println(e => e.Glue(" x ", "C{0}"), "Intersection");
    var tor = set.MaxBy(l => l.Aggregate(1, (acc, i) => acc * i));
    Console.WriteLine($"Ell[{a},{b}](Q) Torsion = {tor!.Glue(" x ", "C{0}")}");
    Console.WriteLine();
}

(int disc, int[]) CandidatsY(int a, int b)
{
    var disc = 4 * a.Pow(3) + 27 * b.Pow(2);
    var r = PrimesDec(int.Abs(disc)).Aggregate(1, (acc, e) => acc * e.Key.Pow(e.Value / 2));
    return (disc, Dividors(r).Append(r).ToArray());
}

IEnumerable<EllPt<Rational>> SolveX(int a, int b, int[] Ys, bool show = false)
{
    var x = FG.BCplxPoly();
    foreach (var y in Ys.Prepend(0))
    {
        var P = x.Pow(3) + a * x + b - y.Pow(2);
        var sols = FG.NRoots(P);
        if (show)
            sols.Println($"Y = {y} P = {P}");

        var ellpts = sols.Where(xi => xi.IsInteger())
            .Select(xi => new EllPt<Rational>($"{(int)double.Round(xi.RealPart.ToDouble)}", $"{y}"));

        foreach (var pt in ellpts)
            yield return pt;
    }
}

void NagellLutz(int a, int b, bool show = false)
{
    var (disc, Ys) = CandidatsY(a, b);
    var ellpts = SolveX(a, b, Ys, show).ToArray();
    var E = new EllGroup<Rational>($"{a}", $"{b}");
    var set = new List<EllPt<Rational>>() { E.O };
    foreach (var pt in ellpts)
    {
        var acc = E.O;
        for (int i = 1; i <= 12; i++)
        {
            acc = E.Op(acc, pt);
            if (acc.IsO)
                break;
        }

        if (acc.IsO)
            set.Add(pt);
    }
    
    if (show)
        set.Println("Elements");
    
    var gEll = Group.Generate(E, set.ToArray());
    DisplayGroup.HeadElements(gEll);
    var abType = Group.AbelianGroupType(gEll);
    Console.WriteLine($"{gEll} Torsion = {abType.Glue(" x ", "C{0}")}");
    Console.WriteLine();
}

// [-36,0]
// https://www.lmfdb.org/EllipticCurve/Q/576/c/3
// 
// [1,0]
// https://www.lmfdb.org/EllipticCurve/Q/64/a/4
// 
// [0,3]
// https://www.lmfdb.org/EllipticCurve/Q/3888/i/2
// 
// [-43,166]
// https://www.lmfdb.org/EllipticCurve/Q/26/b/2
void exampleFp()
{
    Ellmorph(-36, 0, show: true); // Ell[-36,0](Q) Torsion = C2 x C2
    Ellmorph(0, 3); // Ell[0,3](Q) Torsion = C1
    Ellmorph(1, 0); // Ell[1,0](Q) Torsion = C2
    Ellmorph(-43, 166, nbPrimes: 20); // Ell[-43,166](Q) Torsion = C7
}

{
    NagellLutz(-36, 0, show: true); // Ell[-36,0](Q) Torsion = C2 x C2
    NagellLutz(0, 3); // Ell[0,3](Q) Torsion = C1
    NagellLutz(1, 0); // Ell[1,0](Q) Torsion = C2
    NagellLutz(-43, 166); // Ell[-43,166](Q) Torsion = C7
}
// |Ell[-43,166](Q)| = 7
// Type        AbelianGroup
// BaseGroup   Ell[-43,166](Q)
// 
// Elements
// (1)[1] = O
// (2)[7] = (-5,-16)
// (3)[7] = (-5,16)
// (4)[7] = (3,-8)
// (5)[7] = (3,8)
// (6)[7] = (11,-32)
// (7)[7] = (11,32)
// 
// Ell[-43,166](Q) Torsion = C7
// 