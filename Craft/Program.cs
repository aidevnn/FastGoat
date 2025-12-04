using System.Numerics;
using System.Text;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

var seqEq = EqualityComparer<int[]>.Create(
    (a, b) => a is not null && b is not null && a.SequenceEqual(b),
    obj => obj.Length);

IEnumerable<int[]> AllAbTypes(int k)
{
    if (k == 1)
        return [[1]];

    var dec = PrimesDec(k);
    return dec.Select(e => Partitions32[e.Value].Select(l => l.Select(i => e.Key.Pow(i)).ToArray()))
        .MultiLoop()
        .Select(l => l.SelectMany(i => i).OrderDescending().ToArray());
}

(int p, int r)[] AbToElemPows(int[] seq)
{
    return seq.SelectMany(n =>
    {
        if (n < 1)
            throw new();

        if (n == 1)
            return new (int, int)[0];

        return PrimesDec(n).Select(e => (e.Key, e.Value)).Order().ToArray();
    }).Order().ToArray();
}

int[] ElemsToCan(int[] seq)
{
    if (seq.Length == 0 || (seq.Length == 1 && seq[0] == 1))
        return new[] { 1 };

    var set = seq.ToList();
    while (true)
    {
        var arr = set.OrderDescending().ToArray();
        var (e, l) = set.OrderDescending().Select(e => (e, arr.FirstOrDefault(c => Gcd(e, c) == 1, -1)))
            .FirstOrDefault(e => e.Item2 != -1, (-1, -1));
        if (e == -1)
            break;

        set.Remove(e);
        set.Remove(l);
        set.Add(e * l);
    }

    return set.OrderDescending().ToArray();
}

int[] AbToElems(params int[] mods) => AbToElemPows(mods).Select(e => e.p.Pow(e.r)).ToArray();

int[] AbToCan(params int[] mods) => ElemsToCan(AbToElems(mods));

IEnumerable<int[]> Stair(int[] pows)
{
    yield return [0];

    var powsDesc = pows.OrderDescending().ToArray();
    var temp1 = new List<int[]>();
    temp1.AddRange(powsDesc[0].SeqLazy(1).Select(i => new[] { i }).ToArray());
    for (int i = 0; i < powsDesc.Length; i++)
    {
        var temp2 = temp1.ToList();
        temp1.Clear();
        foreach (var current in temp2)
        {
            yield return current;
            if (current.Length == pows.Length)
                continue;

            var max = int.Min(powsDesc[current.Length], current[i]);
            foreach (var j in max.Range(1))
                temp1.Add(current.Append(j).ToArray());
        }
    }
}

IEnumerable<int[]> SubGroups(int[] g)
{
    var g1 = AbToElemPows(g).GroupBy(e => e.p).ToDictionary(
        e => e.Key,
        e => e.Select(f => f.r).OrderDescending().ToArray());

    var gStairs = g1.ToDictionary(e => e.Key, e => Stair(e.Value).ToArray());
    return gStairs.Select(e => e.Value.Select(f => (e.Key, f)).ToArray()).MultiLoop()
        .Select(l =>
            l.Where(li => li.f.All(i => i != 0)).SelectMany(li => li.f.Select(i => li.Key.Pow(i))).ToArray())
        .Select(l => AbToCan(l));
}

IEnumerable<int[]> Morphs(int[] g, int[] n)
{
    var g1Subs = SubGroups(g).ToHashSet(seqEq);
    g1Subs.IntersectWith(SubGroups(n));
    return g1Subs;
}

void TestAbCan(int[] g)
{
    Console.WriteLine($"{g.Glue(" x ", "C{0}")} ~ {AbToCan(g).Glue(" x ", "C{0}")}");
}

void TestAllAbOrd(int o)
{
    var all = AllAbTypes(o).OrderBy(l => l.Length).ThenBy(l => l.Glue(" x ", "C{0:000}")).ToArray();

    all.Println(l => l.Glue(" x ", "C{0}"), $"Abelians groups of order {o} possibles group types :");
    Console.WriteLine($"Total:{all.Length}");
    Console.WriteLine();
}

void TestAllAbOrdDecomp(int o)
{
    var all = AllAbTypes(o).OrderBy(l => l.Length).ThenBy(l => l.Glue(" x ", "C{0:000}")).ToArray();
    var dec = all.Select(g => AbToCan(g)).ToArray();
    var digits = all.Max(g => g.Glue(" x ", "C{0}").Length);
    var fmt = $"{{0,-{digits}}} ~ {{1}}";
    all.Zip(dec).Println(l => string.Format(fmt, l.First.Glue(" x ", "C{0}"), l.Second.Glue(" x ", "C{0}")),
        $"Decomposition of order {o}");
    Console.WriteLine($"Total:{all.Length}");
    Console.WriteLine();
}

void TestStairs(int p, int[] ri)
{
    var th = FG.Abelian(ri.Select(i => p.Pow(i)).ToArray())
        .AllSubgroups().All.Select(e => Group.AbelianGroupType(e).ToArray())
        .DistinctBy(l => l.Glue(" x ", "C{0:000}"))
        .OrderBy(l => l.Length)
        .ThenBy(l => l.Glue(" x ", "C{0:000}"))
        .ToArray();

    var sp = Stair(ri).Select(l => l.Select(i => p.Pow(i)).ToArray())
        .OrderBy(l => l.Length)
        .ThenBy(l => l.Glue(" x ", "C{0:000}"))
        .ToArray();

    var digits = ri.Sum(i => $"C{p.Pow(i)}".Length) + 3 * (ri.Length - 1);
    var fmt = $"{{0,-{digits}}} | {{1,-{digits}}}";

    th.Zip(sp).Println(l => string.Format(fmt, l.First.Glue(" x ", "C{0}"), l.Second.Glue(" x ", "C{0}")));
    Console.WriteLine(th.Length == sp.Length && th.Zip(sp).All(l => l.First.SequenceEqual(l.Second)));
    Console.WriteLine();
}

void TestMorphs(int[] g, int[] n)
{
    var (ag, an) = (FG.Abelian(g), FG.Abelian(n));

    SubGroups(g)
        .OrderBy(l => l.Length)
        .ThenBy(l => l.Glue(" x ", "C{0:000}"))
        .Println(l => l.Glue(" x ", "C{0}"), $"SubGroups(g = {g.Glue(" x ", "C{0}")})");

    SubGroups(n)
        .OrderBy(l => l.Length)
        .ThenBy(l => l.Glue(" x ", "C{0:000}"))
        .Println(l => l.Glue(" x ", "C{0}"), $"SubGroups(n = {n.Glue(" x ", "C{0}")})");

    Group.AllHomomorphisms(ag, an)
        .Select(hom => Group.Generate(an, hom.Image().ToArray()))
        .Select(ab => Group.AbelianGroupType(ab))
        .DistinctBy(l => l, new SequenceEquality<int>())
        .OrderBy(l => l.Length)
        .ThenBy(l => l.Glue(",", "{0:000}"))
        .Println(l => $"g -> {l.Glue(" x ", "C{0}")}", "HomAb(g, n)");

    Morphs(g, n).OrderBy(l => l.Length)
        .ThenBy(l => l.Glue(",", "{0:000}"))
        .Println(l => $"g -> {l.Glue(" x ", "C{0}")}", "HomAb(g, n)");

    Console.WriteLine();
}

void BenchMorphs(int[] g, int[] n)
{
    Action morphsOld = () =>
    {
        var (ag, an) = (FG.Abelian(g), FG.Abelian(n));
        Group.AllHomomorphisms(ag, an)
            .Select(hom => Group.Generate(an, hom.Image().ToArray()))
            .Select(ab => Group.AbelianGroupType(ab))
            .DistinctBy(l => l, new SequenceEquality<int>())
            .ToArray();
    };

    Action morphsNew = () => Morphs(g, n).ToArray();
    
    GlobalStopWatch.Bench(5, "Old", morphsOld);
    GlobalStopWatch.Bench(5, "New", morphsNew);
}

{
    Console.WriteLine("Canonical Decomposition");
    TestAbCan([40, 48]);
    TestAbCan([40, 48, 64]);
    TestAbCan([45, 150, 75]);
    Console.WriteLine();
    
    TestAllAbOrd(450);
    TestAllAbOrd(1600);
    
    TestAllAbOrdDecomp(360);
    TestAllAbOrdDecomp(4320);
    
    TestStairs(2, [1, 2, 2, 4]);
    TestStairs(3, [1, 2, 4]);
    
    TestMorphs([2, 12], [2, 2, 4]);
    TestMorphs([20, 30], [10, 18]);
    
    BenchMorphs([20, 30], [10, 18]);
    BenchMorphs([20, 30], [10, 18]);
    BenchMorphs([20, 30], [10, 18]);
}