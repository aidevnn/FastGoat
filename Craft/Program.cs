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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void TestAbCan(int[] g)
{
    Console.WriteLine($"{g.Glue(" x ", "C{0}")} ~ {AbelianExt.AbToCan(g).Glue(" x ", "C{0}")}");
}

void TestAllAbOrd(int o)
{
    var all = AbelianExt.AllAbTypes(o).OrderBy(l => l.Length).ThenBy(l => l.Glue(" x ", "C{0:000}")).ToArray();

    all.Println(l => l.Glue(" x ", "C{0}"), $"Abelians groups of order {o} possibles group types :");
    Console.WriteLine($"Total:{all.Length}");
    Console.WriteLine();
}

void TestAllAbOrdDecomp(int o)
{
    var all = AbelianExt.AllAbTypes(o).OrderBy(l => l.Length).ThenBy(l => l.Glue(" x ", "C{0:000}")).ToArray();
    var dec = all.Select(g => AbelianExt.AbToCan(g)).ToArray();
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

    var sp = AbelianExt.Stair(ri).Select(l => l.Select(i => p.Pow(i)).ToArray())
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

    AbelianExt.SubGroups(g)
        .OrderBy(l => l.Length)
        .ThenBy(l => l.Glue(" x ", "C{0:000}"))
        .Println(l => l.Glue(" x ", "C{0}"), $"SubGroups(g = {g.Glue(" x ", "C{0}")})");

    AbelianExt.SubGroups(n)
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

    AbelianExt.Morphs(g, n).OrderBy(l => l.Length)
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

    Action morphsNew = () => AbelianExt.Morphs(g, n).ToArray();

    GlobalStopWatch.Bench(5, "Old", morphsOld);
    GlobalStopWatch.Bench(5, "New", morphsNew);
}

void AllTests()
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

{
    // AbelianExt.Quotient([60, 6], [2, 6]).Display();

    int[] n = [20, 30];
    foreach (var s in AbelianExt.SubGroups(n))
        AbelianExt.Quotients(n, s).Display();
    
    Console.WriteLine();

    // var ab = FG.Abelian(n);
    // var subGrs = ab.AllSubgroups();
    // subGrs.Naming();
    // Console.ReadLine();
    // foreach (var sg in subGrs)
    // {
    //     foreach (var g in sg.Conjugates)
    //     {
    //         var q = ab.Over(g);
    //         q.Name = Group.AbelianGroupType(q).Glue(" x ", "C{0}");
    //         DisplayGroup.Head(q);
    //     }
    // }

    // Console.ReadLine();
    // Console.WriteLine(Group.AllHomomorphisms(ab, ab).Count);
}