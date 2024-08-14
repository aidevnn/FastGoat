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
    if (abType.Length <= 3 && g.Count() <= 32)
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

(int, int, int) GenerateAllGroupsOrder(int ord)
{
    if (ord == 96)
        throw new();

    var all = GroupExt.A000001[ord];
    var dividors = Dividors(ord).Where(k => k > 1).ToArray();
    var ab = FG.AllAbelianGroupsOfOrder(ord).Select(e => e.AllSubgroups().ToGroupWrapper());
    var mtCyc = FG.MetaCyclicSdp(ord).Select(e => e.AllSubgroups().ToGroupWrapper());

    var grDiv = dividors.ToDictionary(k => k, k => FG.AllGroupsOfOrder(k).Select(gr => (gr, aut: AutG(gr))));
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

    var lt = ab.AppendIsomorphic(mtCyc, extra, mtCycSdp).Take(all).ToArray();

    var lt1 = lt.Select(e => (subsg: e, names: NamesTree.BuildName(e), rels: Graph.DefiningRelatorsOfGroup(e.Parent)))
        .OrderBy(e => e.subsg.Parent.Count())
        .ThenBy(e => e.subsg.Parent.GroupType)
        .ThenBy(e => e.names[0])
        .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
        .ThenBy(e => e.subsg.Infos)
        .ToArray();

    foreach (var e in lt1)
        Console.WriteLine($"\"{ord};{e.names[0]};{e.rels}\",");

    return (lt.Length, all, all - lt.Length);
}

void AllGoups65to95()
{
    GlobalStopWatch.Restart();
    var l = 31.Range(65).Select(o => (o, m: GenerateAllGroupsOrder(o), t: GroupExt.A000001[o]))
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

    var ord = 96;

    var all = GroupExt.A000001[ord];
    var dividors = Dividors(ord).Where(k => k > 1).ToArray();
    var ab = FG.AllAbelianGroupsOfOrder(ord).Select(e => e.AllSubgroups().ToGroupWrapper());
    var mtCyc = FG.MetaCyclicSdp(ord).Select(e => e.AllSubgroups().ToGroupWrapper());

    var grDiv = dividors.ToDictionary(k => k, k => FG.AllGroupsOfOrder(k).Select(gr => (gr, aut: AutG(gr))));
    var mtCycSdp = dividors.Select(k => (k, m: ord / k)).OrderBy(e => int.Min(e.m, e.k))
        .SelectMany(e =>
            grDiv[e.k].Grid2D(grDiv[e.m]).SelectMany(g => AllSDPFilter(g.t2.gr, g.t1.gr, g.t2.aut, g.t1.aut)))
        .Select(e => e.AllSubgroups().ToGroupWrapper());

    var lt = ab.AppendIsomorphic(mtCyc, mtCycSdp).Take(all).ToArray();
    var founds = lt.SelectMany(e => FG.FindIdGroup(e.Parent, e.Infos).Select(f => f.No)).ToArray();
    var missing = all.Range(1).Except(founds).ToArray();

    FG.AllIds(96).Where(e => missing.Contains(e.No)).Select(e => e.FullName)
        .Println($"Ord{96}:{all} Missing:{missing.Length}");

    lt.Select(e => e.Parent.BaseGroup.Name.Replace("(table)", "")).Distinct().Println();

    GlobalStopWatch.Show(); // ~44m
}

// x: C2
// "Dic12", "GL(2,3)", "D48", "M(6x:8)5", "M(8x:6)3", "M(8x:6)5", "M(8x:6)7", "F(3x:16)2", "M(12x:4)5", "M(12x:4)7", "M(24x:2)5", "M(12x:4)11", "M(24x:2)11", "C6 x D8", "C2 x D24", "C4 x D12", "(C2 x C2) x D12", "C2 x S4", "C4 x A4", "C6 x Q8", "C3 x Q16", "C2 x Dic6", "C2 x SL(2,3)", "(C2 x C2) x A4", "(C2 x C2) x Dic3", "C6 x: D8", "D12 x: C4", "(C4 x C2) x: C6", "(C4 x C4) x: C3", "(C6 x C2) x: C4", "A4 x: C4", "C3 x: Q16", "Dic3 x: C4", "SL(2,3) x: C2", "C3 x: MM16", "(C2 x C2 x C2 x C2) x: C3", "C2 . S4"

// x: C3
// "C8 x C2 x C2", "Q32", "C4 x Q8", "(C2 x C2) x Q8", "(C2 x C2) x: C8", "(C4 x C2) x: C4", "(C8 x C2) x: C2", "C4 x: Q8", "Q8 x: C4", "D8 x: (C2 x C2)", "C8 . C4", "C4 . D8", "(C4 x C4) . C2"

// C3 x:
// "Q32", "M(4x:8)3", "M(8x:4)3", "M(8x:4)5", "M(8x:4)7", "C4 x Q8", "(C4 x C2) x: C4", "C4 x: Q8", "Q8 x: C4", "C8 . C4", "C4 . D8", "(C4 x C4) . C2"

// SL(2,3) x: C4
// (C2 x C2) x: Dic6
// (C2 x C2) x: F(3x:8)2

void NamesGroups96()
{
    GlobalStopWatch.Restart();

    var ord48 = new[]
    {
        "Dic12", "GL(2,3)", "D48", "M(6x:8)5", "M(8x:6)3", "M(8x:6)5", "M(8x:6)7", "F(3x:16)2", "M(12x:4)5",
        "M(12x:4)7", "M(24x:2)5", "M(12x:4)11", "M(24x:2)11", "C6 x D8", "C2 x D24", "C4 x D12", "(C2 x C2) x D12",
        "C2 x S4", "C4 x A4", "C6 x Q8", "C3 x Q16", "C2 x Dic6", "C2 x SL(2,3)", "(C2 x C2) x A4", "(C2 x C2) x Dic3",
        "C6 x: D8", "D12 x: C4", "(C4 x C2) x: C6", "(C4 x C4) x: C3", "(C6 x C2) x: C4", "A4 x: C4", "C3 x: Q16",
        "Dic3 x: C4", "SL(2,3) x: C2", "C3 x: MM16", "(C2 x C2 x C2 x C2) x: C3", "C2 . S4"
    };

    var ord32a = new[]
    {
        "C8 x C2 x C2", "Q32", "C4 x Q8", "(C2 x C2) x Q8", "(C2 x C2) x: C8", "(C4 x C2) x: C4", "(C8 x C2) x: C2",
        "C4 x: Q8", "Q8 x: C4", "D8 x: (C2 x C2)", "C8 . C4", "C4 . D8", "(C4 x C4) . C2"
    };

    var ord32b = new[]
    {
        "Q32", "M(4x:8)3", "M(8x:4)3", "M(8x:4)5", "M(8x:4)7", "C4 x Q8", "(C4 x C2) x: C4", "C4 x: Q8", "Q8 x: C4",
        "C8 . C4", "C4 . D8", "(C4 x C4) . C2"
    };

    var log = Logger.SetOff();
    var ab96 = FG.AllAbelianGroupsOfOrderWg(96).Select(eg => eg.AllSubgroups().ToGroupWrapper()).ToArray();
    var mt96 = FG.MetaCyclicSdpWg(96).Select(eg => eg.AllSubgroups().ToGroupWrapper()).ToArray();
    var gr48 = FG.AllGroupsOfOrder(48).Where(e => ord48.Contains(e.Name)).ToArray();
    var all32 = FG.AllGroupsOfOrder(32).ToArray();
    var gr32a = all32.Where(e => ord32a.Contains(e.Name)).ToArray();
    var gr32b = all32.Where(e => ord32b.Contains(e.Name)).ToArray();
    var sl24 = FG.WordGroup("SL(2,3)", "a4, b3, ababab, a2ba2b-1");
    var dic6 = FG.DiCyclic(6);
    var mt382 = FG.MetaCyclicSdpWg(3, 8, 2);
    var (c2, c3, c4, c2c2) = (FG.AbelianWg(2), FG.AbelianWg(3), FG.AbelianWg(4), FG.AbelianWg(2, 2));
    Logger.Level = log;

    var seq_48_2 = gr48.SelectMany(g => FG.AllSDPFilterLazy(g, c2, true))
        .Select(eg => eg.AllSubgroups().ToGroupWrapper());
    var seq_32_3 = gr32a.SelectMany(g => FG.AllSDPFilterLazy(g, c3, true))
        .Select(eg => eg.AllSubgroups().ToGroupWrapper());
    var seq_3_32 = gr32b.SelectMany(g => FG.AllSDPFilterLazy(c3, g, true))
        .Select(eg => eg.AllSubgroups().ToGroupWrapper());
    var seq_sl23_4 = FG.AllSDPFilterLazy(sl24, c4, true).Select(eg => eg.AllSubgroups().ToGroupWrapper());
    var seq_c2c2_dic6 = FG.AllSDPFilterLazy(c2c2, dic6, true).Select(eg => eg.AllSubgroups().ToGroupWrapper());
    var seq_c2c2_mt382 = FG.AllSDPFilterLazy(c2c2, mt382, true).Select(eg => eg.AllSubgroups().ToGroupWrapper());

    var all = GroupExt.A000001[96];
    var lt0 = ab96.AppendIsomorphic(mt96, seq_48_2, seq_32_3, seq_3_32, seq_sl23_4, seq_c2c2_dic6, seq_c2c2_mt382)
        .Take(all).ToArray();

    var founds = lt0.SelectMany(e => FG.FindIdGroup(e.Parent, e.Infos).Select(f => f.No)).ToArray();
    var missing = all.Range(1).Except(founds).ToArray();

    FG.AllIds(96).Where(e => missing.Contains(e.No)).Select(e => e.FullName)
        .Println($"Ord{96}:{all} Missing:{missing.Length}");

    GlobalStopWatch.Show();

    // Group (C2 x C2 x D12) x: C2 has 441 normals subgroups
    var extra = FG.WordGroup("(C2 x C2) x D12", "a6, b2, c2, d2, abab, bcbc, bdbd, cdcd, caca-1, dada-1")
        .ToGroupWrapper();
    var xtrLeaf = new Leaf(extra, extra.Name);
    var c2Leaf = new Leaf(c2.ToGroupWrapper());
    log = Logger.SetOff();
    var lt1 = lt0.Select((e, k) =>
        {
            Console.Write($"no{k + 1:000}, ");
            var rels = Graph.DefiningRelatorsOfGroup(e.Parent);
            if (e.Parent.BaseGroup.Name.Contains(extra.Name) && e.Infos.AllNorms == 441)
            {
                var name0 = new SemiDirectProductOp(xtrLeaf, c2Leaf, e.Parent);
                return (subsg: e, names: [name0], rels);
            }
            else
                return (subsg: e, names: NamesTree.BuildName(e), rels);
        })
        .OrderBy(e => e.subsg.Parent.Count())
        .ThenBy(e => e.subsg.Parent.GroupType)
        .ThenBy(e => e.names[0])
        .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
        .ThenBy(e => e.subsg.Infos)
        .ToArray();

    Logger.Level = log;
    Console.WriteLine();
    
    foreach (var e in lt1)
        Console.WriteLine($"\"{96};{e.names[0]};{e.rels}\",");

    Console.WriteLine();
    GlobalStopWatch.Show("End");
}

void Order81()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();
    var ab81 = FG.AllAbelianGroupsOfOrder(81).Select(e => e.AllSubgroups().ToGroupWrapper());
    var ext81 = FG.AllExtensions((FG.Abelian(9, 3), FG.Abelian(3))).Select(e => e.allSubs.ToGroupWrapper());
    var sdp81 = FG.AllSDPFilterLazy(FG.ElementaryAbelian(27), FG.Abelian(3))
        .Select(e => e.AllSubgroups().ToGroupWrapper());
    ab81.AppendIsomorphic(ext81, sdp81)
        .Take(GroupExt.A000001[81])
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
}

{
    // AllGoups65to95();

    Logger.Level = LogLevel.Level1;
    Group.ActivedStorage(false);

    // AllGroups96(); // Time:41m26s
    // NamesGroups96(); // Time:27m26s
    
    var sum0 = GroupExt.A000001.Take(97).Sum();
    var sum1 = GroupExt.DB.Length;
    Console.WriteLine(new { sum0, sum1 }); // { sum0 = 1024, sum1 = 1024 }
}