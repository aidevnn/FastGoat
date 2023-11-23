using System.Collections;
using System.ComponentModel;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
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
var nbOps = 10000;

IEnumerable<(CrMap<Tn, Tg> c, ConcreteGroup<Ep2<Tn, Tg>> ext, Dictionary<ConcreteGroup<Ep2<Tn, Tg>>,
    List<ConcreteGroup<Ep2<Tn, Tg>>>> allSubs, (int, int, int) infos)> AllExtensions<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var CN = Group.Zentrum(N);
    var autN = Group.AutomorphismGroup(N);
    var ops = Group.AllHomomorphisms(G, autN);
    var dicExts = new Dictionary<(int, int, int), HashSet<ConcreteGroup<Ep2<Tn, Tg>>>>();
    foreach (var (op, i) in ops.Select((l, i) => (l, i + 1)).Take(nbOps))
    {
        var L = op.ToMapElt(autN);
        var lbl = $"Lbl{i}/{ops.Count}";
        var (cohs, cobs, cocs) = ZNSolver.ReduceCohomologies(CN, G, L, lbl: lbl);
        foreach (var c in cohs)
        {
            var c0 = c.ToMapElt;
            var ext = Group.ExtensionGroup(N, L, c0, G);
            var extBase = ((ExtensionGroupBase<Tn, Tg>)ext.BaseGroup);
            if (extBase.IsGroup)
            {
                var allSubs = Group.AllSubGroups(ext);
                var infos = CocyclesDFS.SubGroupsDetails(allSubs);
                if (dicExts.ContainsKey(infos))
                {
                    if (dicExts[infos].Add(ext))
                    {
                        NameGroup(allSubs);
                        var g = allSubs.MaxBy(e => e.Key.Count()).Key;
                        ext.SetName(g.Name);
                        yield return (c, ext, allSubs, infos);
                    }
                }
                else
                {
                    dicExts[infos] = new(new IsomorphEquality<Ep2<Tn, Tg>>()) { ext };
                    NameGroup(allSubs);
                    var g = allSubs.MaxBy(e => e.Key.Count()).Key;
                    ext.SetName(g.Name);
                    yield return (c, ext, allSubs, infos);
                }
            }
            else
            {
                Console.WriteLine("????????????????????? Extension isnt a group");
            }
        }

        Console.WriteLine($"Nb Exts:{dicExts.Values.Sum(v => v.Count)}");
    }
}


IEnumerable<(CrMap<Tn, Tg> c, ConcreteGroup<Ep2<Tn, Tg>> ext, Dictionary<ConcreteGroup<Ep2<Tn, Tg>>,
        List<ConcreteGroup<Ep2<Tn, Tg>>>> allSubs, (int, int, int) infos)>
    AllExtensions2<Tn, Tg>(params (ConcreteGroup<Tn>, ConcreteGroup<Tg>)[] tuples)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var dicExts = new Dictionary<(int, int, int), HashSet<ConcreteGroup<Ep2<Tn, Tg>>>>();
    foreach (var (n, g) in tuples)
    {
        foreach (var extInfos in AllExtensions(n, g))
        {
            if (dicExts.ContainsKey(extInfos.infos))
            {
                if (dicExts[extInfos.infos].Add(extInfos.ext))
                    yield return extInfos;
            }
            else
            {
                dicExts[extInfos.infos] = new(new IsomorphEquality<Ep2<Tn, Tg>>()) { extInfos.ext };
                yield return extInfos;
            }
        }

        Console.WriteLine();
        Console.WriteLine($"Total Exts:{dicExts.Values.Sum(v => v.Count)}");
    }
}

Dictionary<ConcreteGroup<T>, ConcreteGroup<T>[]>
    IsProduct<T>(Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> subGroups)
    where T : struct, IElt<T>
{
    var g = subGroups.MaxBy(e => e.Key.Count()).Key;
    var og = g.Count();
    var tr = Group.Generate("Tr", g, g.Neutral());
    var normalSubGroups = subGroups.Where(e => e.Value.Count == 1).Select(e => e.Key).ToArray();
    var otherSubroups = subGroups.OrderBy(e => e.Value.Count).SelectMany(e => e.Value).Except(normalSubGroups).ToArray();
    var properNormalSubGroups = normalSubGroups.Where(e => e.Count() != 1 && e.Count() != og)
        .OrderBy(e => e.GroupType == GroupType.NonAbelianGroup)
        .ThenByDescending(e => e.Count())
        .ToArray();
    return properNormalSubGroups.ToDictionary(n => n,
            n => otherSubroups.OrderByDescending(e => subGroups.Any(kv => kv.Key.SetEquals(e)) ? e.Count() : 0)
                .ThenBy(e => e.Count())
                .Where(e => e.Count() * n.Count() == og && n.Intersect(e).Count() == 1).ToArray())
        .Where(e => e.Value.Length != 0).ToDictionary(e => e.Key, e => e.Value);
}

Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>>
    SubGroupRes<T>(ConcreteGroup<T> g, Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> subGroups)
    where T : struct, IElt<T>
{
    return subGroups.Values.Select(e => e.Where(e0 => e0.SubSetOf(g)).ToList()).Where(e => e.Count > 0)
        .ToDictionary(e => e.First(), e => e);
}

void NameGroup<T>(Dictionary<ConcreteGroup<T>, List<ConcreteGroup<T>>> subGroups)
    where T : struct, IElt<T>
{
    var g = subGroups.MaxBy(e => e.Key.Count()).Key;
    if (g.GroupType == GroupType.AbelianGroup)
    {
        var name = Group.AbelianInvariants(g).Select(e => e.o).Glue(" x ", "C{0}");
        g.SetName(name);
    }
    else
    {
        var dic = IsProduct(subGroups);
        if (dic.Count != 0)
        {
            var (k, h) = dic.Select(e => (e.Key, e.Value[0]))
                .OrderByDescending(e => subGroups.Any(kv => kv.Key.SetEquals(e.Key)) &&
                                        subGroups.Any(kv => kv.Key.SetEquals(e.Item2))) // direct product
                .ThenByDescending(e => e.Key.GroupType == GroupType.AbelianGroup && e.Item2.Count() == 2) // dihedral
                .ThenByDescending(e => e.Key.GroupType == GroupType.AbelianGroup &&
                                       e.Key.Name.Count(c => c == 'C') == 1 &&
                                       e.Item2.GroupType == GroupType.AbelianGroup &&
                                       e.Item2.Name.Count(c => c == 'C') == 1) // pq-group
                .ThenByDescending(e => subGroups.Any(kv => kv.Key.SetEquals(e.Item2)) ? e.Item2.Count() : e.Key.Count())
                .First();
            var subGroupsK = SubGroupRes(k, subGroups);
            var subGroupsH = SubGroupRes(h, subGroups);

            NameGroup(subGroupsK);
            var gk = subGroupsK.MaxBy(e => e.Key.Count()).Key;
            var nk = gk.Name.Contains(' ') ? $"({gk.Name})" : gk.Name;
            NameGroup(subGroupsH);
            var gh = subGroupsH.MaxBy(e => e.Key.Count()).Key;
            var nh = gh.Name.Contains(' ') ? $"({gh.Name})" : gh.Name;
            var sep = subGroups.Any(e => e.Value.Count == 1 && e.Key.SetEquals(h)) ? "x" : "x:";
            g.SetName($"{nk} {sep} {nh}");
        }
        else
        {
            g.SetName($"G[{g.Count()}]");
        }
    }
}

void Order16()
{
    GlobalStopWatch.Restart();
    nbOps = 4;
    var allExts16 = AllExtensions2(TestTwoCohomology.AllAbelianGroupsOrder(8).Select(e => (e, FG.Abelian(2))).ToArray())
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    CocyclesDFS.DisplayInfosGroups(allExts16.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext16:{allExts16.Length}");
    Console.Beep();
}

void Order24()
{
    GlobalStopWatch.Restart();
    nbOps = 3;
    var tuples8 = TestTwoCohomology.AllAbelianGroupsOrder(4).Select(e => (e, FG.Abelian(2)));
    var allExts8 = AllExtensions2(tuples8.ToArray()).ToArray();

    var tuples12a = TestTwoCohomology.AllAbelianGroupsOrder(4).Select(e => (e, FG.Abelian(3)));
    var tuples12b = TestTwoCohomology.AllAbelianGroupsOrder(6).Select(e => (e, FG.Abelian(2)));
    var tuples12 = tuples12a.Concat(tuples12b).ToArray();
    var allExts12 = AllExtensions2(tuples12).ToArray();

    var tuples24a = allExts8.Select(e => (e.ext, FG.Abelian(3)));
    var tuples24b = allExts12.Select(e => (e.ext, FG.Abelian(2)));
    var tuples24 = tuples24a.Concat(tuples24b).ToArray();
    var allExts24 = AllExtensions2(tuples24)
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    CocyclesDFS.DisplayInfosGroups(allExts8.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext8:{allExts8.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts12.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext12:{allExts12.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts24.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext24:{allExts24.Length}");
    Console.Beep();
}

void Order42()
{
    GlobalStopWatch.Restart();
    nbOps = 4;
    var tuples14 = TestTwoCohomology.AllAbelianGroupsOrder(7).Select(e => (e, FG.Abelian(2)));
    var allExts14 = AllExtensions2(tuples14.ToArray()).ToArray();
    var tuples21 = TestTwoCohomology.AllAbelianGroupsOrder(7).Select(e => (e, FG.Abelian(3)));
    var allExts21 = AllExtensions2(tuples21.ToArray()).ToArray();

    var tuples42a = allExts14.Select(e => (e.ext, FG.Abelian(3)));
    var tuples42b = allExts21.Select(e => (e.ext, FG.Abelian(2)));
    var tuples42 = tuples42a.Concat(tuples42b).ToArray();
    var allExts42 = AllExtensions2(tuples42)
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    CocyclesDFS.DisplayInfosGroups(allExts14.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext14:{allExts14.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts21.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext21:{allExts21.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts42.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext42:{allExts42.Length}");
    Console.Beep();
}

void Order40()
{
    GlobalStopWatch.Restart();
    nbOps = 4;
    var allExts10 = AllExtensions2((FG.Abelian(5), FG.Abelian(2))).ToArray();
    var tuples20 = TestTwoCohomology.AllAbelianGroupsOrder(4).Select(e => (FG.Abelian(5), e)).ToArray();
    var allExts20 = AllExtensions2(tuples20).ToArray();
    var tuples40 = TestTwoCohomology.AllAbelianGroupsOrder(4).Select(e => (FG.Abelian(2, 5), e)).ToArray();
    var allExts40 = AllExtensions2(tuples40)
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    CocyclesDFS.DisplayInfosGroups(allExts10.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext10:{allExts10.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts20.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext20:{allExts20.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts40.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext40:{allExts40.Length}");
    Console.Beep();
}

void Order32()
{
    GlobalStopWatch.Restart();
    var exts = AllExtensions2(
            (FG.Abelian(16), FG.Abelian(2)),
            (FG.Abelian(2, 8), FG.Abelian(2)),
            (FG.Abelian(4, 4), FG.Abelian(2)),
            (FG.Abelian(2, 2, 4), FG.Abelian(2)),
            (FG.Abelian(2, 4), FG.Abelian(4)),
            (FG.Abelian(2, 2, 2), FG.Abelian(4)),
            (FG.Abelian(8), FG.Abelian(2, 2)),
            (FG.Abelian(2, 4), FG.Abelian(2, 2)),
            (FG.Abelian(2, 2, 2), FG.Abelian(2, 2)))
        .Take(51)
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    CocyclesDFS.DisplayInfosGroups(exts.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext:{exts.Length}"); // # Nb Ext:51 Time:781578 ms ~ 13min

    Console.Beep();
}

void Order56()
{
    GlobalStopWatch.Restart();
    nbOps = 2;
    var allExts28 = AllExtensions2((FG.Abelian(7), FG.Abelian(4)), (FG.Abelian(7), FG.Abelian(2, 2)))
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    var allExts56 = AllExtensions2(
            (FG.Abelian(2, 7), FG.Abelian(4)),
            (FG.Abelian(2, 7), FG.Abelian(2, 2)),
            (FG.Abelian(2, 2, 2), FG.Abelian(7)))
        .Take(13)
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    CocyclesDFS.DisplayInfosGroups(allExts28.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext56:{allExts28.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts56.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext56:{allExts56.Length}");
    Console.Beep();
}

void Order81()
{
    GlobalStopWatch.Restart();
    nbOps = 2;
    var allExts27 = AllExtensions2((FG.Abelian(9), FG.Abelian(3)), (FG.Abelian(3, 3), FG.Abelian(3)))
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    (allExts27[2], allExts27[4]) = (allExts27[4], allExts27[2]);
    var allExts81 = AllExtensions2(allExts27.Select(e => (e.ext, FG.Abelian(3))).ToArray())
        .Take(15)
        .OrderBy(e => e.ext.GroupType)
        .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
        .ThenBy(e => e.infos).ToArray();

    CocyclesDFS.DisplayInfosGroups(allExts27.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext27:{allExts27.Length}");
    CocyclesDFS.DisplayInfosGroups(allExts81.Select(e => (e.ext, e.infos)).ToArray(), naming: false);
    GlobalStopWatch.Show($"Nb Ext81:{allExts81.Length}"); // # Nb Ext81:15 Time:489456 ms

    Console.Beep();
}

{
    // Order16();
    // Order24();
    // Order32();
    // Order40();
    // Order42();
    Order56();
    // Order81();
}