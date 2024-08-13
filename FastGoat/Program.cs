using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
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

ConcreteGroup<Automorphism<T>> AutG<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    if (g.GroupType == GroupType.NonAbelianGroup)
        return Group.AutomorphismGroup(g);

    var abType = Group.AbelianGroupType(g);
    if (abType.Length <= 3 || (abType.Length == 4 && abType.All(k => k == 2)))
        return Group.AutomorphismGroup(g);

    var autG = Group.AutBase(g);
    return Group.Generate($"Aut[{g}]", autG, autG.Neutral());
}

IEnumerable<SemiDirectProduct<T1, T2>> AllSDPFilter<T1, T2>(ConcreteGroup<T1> N, ConcreteGroup<T2> G,
    ConcreteGroup<Automorphism<T1>> autN, ConcreteGroup<Automorphism<T2>> autG)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    if (Logger.Level != LogLevel.Off)
        Console.WriteLine($"############### AllSDP {N.NameParenthesis()} x: {G.NameParenthesis()}");

    var allOps = Group.AllHomomorphisms(G, autN);
    if (Logger.Level != LogLevel.Off)
        Console.WriteLine($"AutG:{autG.Count()} AutN:{autN.Count()}");

    var ops = allOps.Distinct(FG.EqOpByAut(G, autG, autN));
    var k = 1;
    foreach (var theta in ops)
    {
        if (Logger.Level != LogLevel.Off)
            Console.WriteLine($"  ##   {k++,3}   ##");

        yield return Group.SemiDirectProd(N, theta, G);
    }
}

(int, int, int) AbelianMetacyclics(int ord)
{
    var all = GroupExt.A000001[ord];
    if (ord == 81) --all;
    var dividors = Dividors(ord).Where(k => k > 1).ToArray();
    var ab = FG.AllAbelianGroupsOfOrder(ord).Select(e => e.AllSubgroups().ToGroupWrapper());
    var mtCyc = FG.MetaCyclicSdp(ord).Select(e => e.AllSubgroups().ToGroupWrapper());

    var grDiv = dividors.ToDictionary(k => k, k => FG.AllGroupsOfOrder(k).Select(gr => (gr, aut: AutG(gr))));
    var mtCycProd = grDiv.SelectMany(e =>
            e.Value.Grid2D(grDiv[ord / e.Key]).Select(g => Product.Generate(g.t1.gr, g.t2.gr)))
        .Select(e => e.AllSubgroups().ToGroupWrapper());
    var mtCycSdp = dividors.Select(k => (k, m: ord / k)).OrderBy(e => int.Min(e.m, e.k))
        .SelectMany(e =>
            grDiv[e.k].Grid2D(grDiv[e.m]).SelectMany(g => AllSDPFilter(g.t2.gr, g.t1.gr, g.t2.aut, g.t1.aut)))
        .Select(e => e.AllSubgroups().ToGroupWrapper());
    // var mtCycSdp = grDiv.SelectMany(e =>
    //         e.Value.Grid2D(grDiv[ord / e.Key])
    //             .SelectMany(g => AllSDPFilter(g.t2.gr, g.t1.gr, g.t2.aut, g.t1.aut)))
    //     .Select(e => e.AllSubgroups().ToGroupWrapper());

    Console.WriteLine("##### Go");
    var lt = ab.AppendIsomorphic(mtCyc, mtCycSdp).Take(all).ToArray();
    var founds = lt.SelectMany(e => FG.FindIdGroup(e.Parent, e.Infos).Select(f => f.No)).ToArray();
    var missing = all.Range(1).Except(founds).ToArray();

    FG.AllIds(ord).Where(e => missing.Contains(e.No)).Select(e => e.FullName)
        .Println($"Ord{ord}:{all} Missing:{missing.Length}");

    return (lt.Length, all, missing.Length);
}

void AllGoups65to95()
{
    GlobalStopWatch.Restart();
    var l = 31.Range(65).Select(o => (o, m: AbelianMetacyclics(o), t: GroupExt.A000001[o]))
        .ToArray();

    l.Where(k => k.m.Item3 != 0).Println();
    Console.WriteLine(l.Select(e => e.m.Item1).Sum());
    Console.WriteLine(l.Select(e => e.m.Item2).Sum());
    Console.WriteLine(l.Select(e => e.m.Item3).Sum());
    GlobalStopWatch.Show(); // 1m40s
}

void AllGroups96()
{
    GlobalStopWatch.Restart();
    Group.ActivedStorage(false);
    AbelianMetacyclics(96);
    GlobalStopWatch.Show(); // ~44m
}

{
    Logger.Level = LogLevel.Level1;
    AllGroups96();
    // AllGoups65to95();
}