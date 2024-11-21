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
    Console.WriteLine($"E[{a},{b}](Q) Torsion = {tor!.Glue(" x ", "C{0}")}");
    Console.WriteLine();
}

{
    Ellmorph(-36, 0, show: true); // E[-36,0](Q) Torsion = C2 x C2
    Ellmorph(0, 3); // E[0,3](Q) Torsion = C1
    Ellmorph(1, 0); // E[1,0](Q) Torsion = C2
    Ellmorph(-43, 166, nbPrimes: 20); // E[-43,166](Q) Torsion = C7
}