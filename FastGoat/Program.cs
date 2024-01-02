using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

{
    GlobalStopWatch.Restart();
    
    foreach (var g in BatchGroups())
        TestAllSubgroups(g);
    
    foreach (var g in BatchGroups())
        BenchAllSubgroups(g);
    
    TestAllSubgroups(FG.ElementaryAbelian(64));
}

IEnumerable<ConcreteGroup<TableElt>> BatchGroups()
{
    for (int k = 3; k < 10; k++)
        yield return FG.Dihedral(k).ToCGTable();

    for (int k = 3; k < 10; k++)
        yield return FG.DiCyclic(k).ToCGTable();

    for (int k = 3; k < 7; k++)
        yield return FG.Quaternion(2.Pow(k)).ToCGTable();

    for (int o = 2; o <= 32; o++)
    {
        foreach (var ab in FG.AllAbelianGroupsOfOrder(o))
            yield return ab.ToCGTable();
    }
    
    for (int n = 3; n <= 6; n++)
    {
        yield return FG.Alternate(n).ToCGTable();
        yield return FG.Symmetric(n).ToCGTable();
    }
}

void BenchAllSubgroups<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    Console.WriteLine(g);
    GlobalStopWatch.Bench(5, "New Meth", () => AllSubGroupsNewMeth(g));
    GlobalStopWatch.Bench(5, "Old Meth", () => g.AllSubgroups());
    Console.WriteLine();
}

void TestAllSubgroups<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    Console.WriteLine(g);
    GlobalStopWatch.AddLap();
    var subs = AllSubGroupsNewMeth(g);
    var all = subs.Sum(e => e.Value.Count);
    var norms = subs.Count(e => e.Value.Count == 1);
    var actual = new SubGroupsInfos(all, subs.Count, norms);
    GlobalStopWatch.Show($"actual:  {actual}");
    GlobalStopWatch.AddLap();
    var expected = g.AllSubgroups().Infos;
    GlobalStopWatch.Show($"expected:{expected}");
    Console.WriteLine();
    
    if (actual != expected)
        throw new();
}

Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> AllSubGroupsNewMeth<T>(ConcreteGroup<T> g)
    where T : struct, IElt<T>
{
    var og = g.Count();
    var allSubGrs = new HashSet<ConcreteGroup<T>>(og * og, new GroupSetEquality<T>());
    var table = new Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>>(og * og, new GroupSetEquality<T>());

    foreach (var e0 in g)
    {
        var cyc = Group.Generate(g, e0);
        if (!allSubGrs.Contains(cyc))
        {
            var conjs = table[cyc] = Group.SubGroupsConjugates(g, cyc);
            allSubGrs.UnionWith(conjs);
        }
    }

    var sgsRem = table.Keys.ToHashSet(new GroupSetEquality<T>());
    var cycRem = allSubGrs.ToHashSet(new GroupSetEquality<T>());
    while (sgsRem.Count != 0)
    {
        var sgs = sgsRem.OrderBy(sg => sg.Count()).ToArray();
        var cyc = cycRem.OrderBy(sg => sg.Count()).ToArray();
        sgsRem.Clear();
        cycRem.Clear();
        foreach (var (sg0, sg1) in cyc.Grid2D(sgs))
        {
            if (sg0.SuperSetOf(sg1) || sg0.SubSetOf(sg1))
                continue;

            var gens = sg0.GetGenerators().Union(sg1.GetGenerators()).ToArray();
            var elts = Group.GenerateElements(g, sg1.Concat(sg0).ToHashSet(), gens.ToList()).ToHashSet();
            if (allSubGrs.All(g0 => !g0.SetEquals(elts)))
            {
                var sg2 = Group.Generate(g, gens);
                var conjsSg2 = table[sg2] = Group.SubGroupsConjugates(g, sg2);
                allSubGrs.UnionWith(conjsSg2);
                sgsRem.Add(sg2);
                cycRem.Add(sg0);
            }
        }
    }

    var gName = g.NameParenthesis();
    foreach (var (g0, i) in table.Keys.OrderBy(g0 => g0.Count()).Select((g0, i) => (g0, i)))
    {
        g0.Name = $"{gName}-SubGr{i + 1}";
        foreach (var (cg0, j) in table[g0].Select((cg0, j) => (cg0, j)))
            cg0.Name = $"{g0}-Cj{j + 1}";
    }

    var (k, v) = table.First(g0 => g0.Key.Count() == og);
    k.Name = g.Name;
    v.First().Name = g.Name;
    return table;
}