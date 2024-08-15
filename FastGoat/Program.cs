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
    if (abType.Length <= 3)
        return Group.AutomorphismGroup(g);
    else
    {
        var autG = Group.AutBase(g);
        return Group.Generate($"Aut[{g}]", autG, autG.Neutral());
    }
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

AllSubgroups<WElt>[] GenerateAllGroupsOrder(int ord)
{
    var all = GroupExt.A000001[ord];
    var dividors = Dividors(ord).Where(k => k > 1).ToArray();
    var ab = FG.AllAbelianGroupsOfOrder(ord).Select(e => e.AllSubgroups().ToGroupWrapper());
    var mtCyc = FG.MetaCyclicSdp(ord).Select(e => e.AllSubgroups().ToGroupWrapper());

    var grDiv = dividors.ToDictionary(
        k => k,
        k => FG.AllGroupsOfOrder(k).Where(g => !g.Name.Contains("x C3 x C3")).Select(gr => (gr, aut: AutG(gr))));
    var mtCycSdp = dividors.Select(k => (k, m: ord / k)).OrderBy(e => int.Min(e.m, e.k))
        .SelectMany(e =>
            grDiv[e.k].Grid2D(grDiv[e.m]).SelectMany(g => AllSDPFilter(g.t2.gr, g.t1.gr, g.t2.aut, g.t1.aut)))
        .Select(e => e.AllSubgroups().ToGroupWrapper());

    var extra = new List<AllSubgroups<WElt>>();
    if (ord == 80)
        extra.Add(FG.WordGroup("a5,b2,c2,d2,e2,bcbc,bdbd,bebe,cdcd,cece,dede,ba-1ca,ca-1da,da-1ea,bceaeca-1")
            .AllSubgroups().ToGroupWrapper());
    if (ord == 81)
        extra.Add(FG.WordGroup("c3,a3b3,aca2c-1,cbc-1b-1,abca-1b-1").AllSubgroups().ToGroupWrapper());
    if (ord == 120)
        extra.Add(FG.WordGroup("SL(2,5)", "a4,b3,a2ba2b-1,ababababa-1b").AllSubgroups().ToGroupWrapper());

    var lt = ab.AppendIsomorphic(mtCyc, extra, mtCycSdp).Take(all).ToArray();
    var founds = lt.SelectMany(e => FG.FindIdGroup(e.Parent, e.Infos).Select(f => f.No)).ToArray();
    var missing = all.Range(1).Except(founds).ToArray();

    FG.AllIds(ord).Where(e => missing.Contains(e.No)).Select(e => e.FullName)
        .Println($"Ord{ord}:{all} Missing:{missing.Length}");

    Console.WriteLine();
    return lt;
}

void AllGoupsStart(int start, int nb = 1)
{
    GlobalStopWatch.Restart();
    var l = nb.Range(start).Select(o => (o, m: GenerateAllGroupsOrder(o), t: GroupExt.A000001[o])).ToArray();

    l.Where(k => k.m.Length != 0).Println(e => $"ord:{e.o} all:{e.t} found:{e.m.Length} missing:{e.t - e.m.Length}", "Summary");
    Console.WriteLine(l.Select(e => e.m.Length).Sum());
    Console.WriteLine(l.Select(e => e.t).Sum());
    Console.WriteLine(l.Select(e => e.t - e.m.Length).Sum());

    GlobalStopWatch.Show();
    Console.WriteLine();

    var log = Logger.SetOff();
    foreach (var (ord, lt, _) in l)
    {
        var lt1 = lt.Select(e =>
                (subsg: e, names: NamesTree.BuildName(e), rels: Graph.DefiningRelatorsOfGroup(e.Parent)))
            .OrderBy(e => e.subsg.Parent.Count())
            .ThenBy(e => e.subsg.Parent.GroupType)
            .ThenBy(e => e.names[0])
            .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
            .ThenBy(e => e.subsg.Infos)
            .ToArray();

        foreach (var e in lt1)
            Console.WriteLine($"\"{ord};{e.names[0]};{e.rels}\",");
    }

    Logger.Level = log;

    GlobalStopWatch.Show();
}

{
    Logger.Level = LogLevel.Level1;
    
    AllGoupsStart(65, 31); // Time:1m44s
    
    // AllGoupsStart(96, 1); // Time:23m52s
    
    // AllGoupsStart(97, 31); // Time:4m45s
    
    Console.Beep();
}
