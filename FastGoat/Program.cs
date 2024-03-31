using System.Collections;
using System.ComponentModel;
using System.ComponentModel.DataAnnotations;
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
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void GroupsFromDB(int max = 16)
{
    GlobalStopWatch.Restart();

    GroupExt.DB.Select(s => s.Split(';'))
        .Where(s => int.Parse(s[0]) <= max)
        .Select(s =>
        {
            Logger.Level = LogLevel.Off;
            var g = FG.WordGroup(s[1], s[2]);
            Logger.Level = LogLevel.Level1;
            return g;
        })
        .Select(g => g.AllSubgroups().ToGroupWrapper())
        .FilterIsomorphic()
        .Naming()
        .DisplayNames(showBasegroup: false)
        .CheckMissings();

    GlobalStopWatch.Show("End");
}

IEnumerable<AllSubgroups<WElt>> GetAllProds(int ord)
{
    var allAb = FG.AllAbelianGroupsOfOrder(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    var allMetaCyc = FG.MetaCyclicSdp(ord).Select(g => g.AllSubgroups().ToGroupWrapper());
    return allAb.Concat(allMetaCyc);
}

(AllSubgroups<WElt> subsg, ANameElt[] names)[] AllGroupsUpTo32()
{
    var allAb = 32.Range(1).SelectMany(k => FG.AllAbelianGroupsOfOrder(k)).Select(e => e.AllSubgroups().ToGroupWrapper());
    var allMtCycSdp = 32.Range(1).SelectMany(k => FG.MetaCyclicSdp(k)).Select(e => e.AllSubgroups().ToGroupWrapper());
    var ext8 = FG.AllExtensions((FG.Abelian(2), FG.Abelian(2, 2))).Select(e => e.allSubs.ToGroupWrapper());
    var ext12 = FG.AllExtensions((FG.Abelian(2, 2), FG.Abelian(3))).Select(e => e.allSubs.ToGroupWrapper());
    var ext16 = FG.AllExtensions((FG.Abelian(8), FG.Abelian(2)), (FG.Abelian(4, 2), FG.Abelian(2)))
        .Select(e => e.allSubs.ToGroupWrapper());
    var ext18 = FG.AllExtensions((FG.Abelian(3, 3), FG.Abelian(2))).Select(e => e.allSubs.ToGroupWrapper());
    var ext24 = FG.AllExtensions(
            (FG.Abelian(12).ToCGW(), FG.Abelian(2)),
            (FG.Abelian(2, 6).ToCGW(), FG.Abelian(2)),
            (FG.Alternate(4).ToCGW(), FG.Abelian(2)),
            (FG.Quaternion(8).ToCGW(), FG.Abelian(3)))
        .Select(e => e.allSubs.ToGroupWrapper());
    var ext27 = FG.AllExtensions((FG.Abelian(3, 3), FG.Abelian(3))).Select(e => e.allSubs.ToGroupWrapper());
    var ext32 = FG.AllExtensions(
            (FG.Abelian(16), FG.Abelian(2)),
            (FG.Abelian(8, 2), FG.Abelian(2)),
            (FG.Abelian(4, 4), FG.Abelian(2)),
            (FG.Abelian(4, 2), FG.Abelian(4)),
            (FG.Abelian(4, 2), FG.Abelian(2, 2)))
        .Select(e => e.allSubs.ToGroupWrapper());

    return allAb.AppendIsomorphic(allMtCycSdp, ext8, ext12, ext16, ext18, ext27, ext24, ext32)
        .Naming()
        .ToArray();
}

IEnumerable<AllSubgroups<WElt>> Ord48()
{
    var ord24 = FG.AllExtensions(
            (FG.Abelian(12).ToCGW(), FG.Abelian(2)),
            (FG.Abelian(2, 6).ToCGW(), FG.Abelian(2)),
            (FG.Alternate(4).ToCGW(), FG.Abelian(2)),
            (FG.Quaternion(8).ToCGW(), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .FilterIsomorphic()
        .Take(GroupExt.A000001[24])
        .Naming()
        .ToArray();

    var listByC2 = ord24.SelectMany(e => FG.AllExtensions((e.subsg.Parent, FG.Abelian(2)))).Select(e => e.allSubs.ToGroupWrapper());
    var listByC3 = FG.AllAbelianGroupsOfOrder(16)
        .SelectMany(e => FG.AllSDPFilterLazy(e, FG.Abelian(3)))
        .Select(e => e.AllSubgroups().ToGroupWrapper());

    var ord = 48;
    return listByC2.AppendIsomorphic(listByC3, GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord36()
{
    var ord = 36;
    return FG.AllExtensions(
            (FG.Abelian(9), FG.Abelian(4)),
            (FG.Abelian(9), FG.Abelian(2, 2)),
            (FG.Abelian(3, 3), FG.Abelian(4)),
            (FG.Abelian(3, 3), FG.Abelian(2, 2)),
            (FG.Abelian(2, 6), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord40()
{
    var ord = 48;
    return FG.AllExtensions(
            (FG.Abelian(10), FG.Abelian(4)),
            (FG.Abelian(10), FG.Abelian(2, 2))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord50()
{
    var ord = 50;
    return FG.AllSDPFilter(FG.Abelian(5, 5), FG.Abelian(2))
        .Select(g => g.AllSubgroups().ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord54()
{
    var ord = 54;
    return FG.AllExtensions(
            (FG.Abelian(3, 6), FG.AbelianPerm(3)),
            (FG.Abelian(3, 3), FG.Symmetric(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord56()
{
    var ord = 56;
    var sdp = FG.AllSDPFilter(FG.Abelian(2, 2, 2), FG.Abelian(7)).Select(g => g.AllSubgroups().ToGroupWrapper());
    return FG.AllExtensions(
            (FG.Abelian(28), FG.Abelian(2)),
            (FG.Abelian(2, 14), FG.Abelian(2))
        )
        .Select(e => e.allSubs.ToGroupWrapper())
        .AppendIsomorphic(sdp, GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

IEnumerable<AllSubgroups<WElt>> Ord60()
{
    var ord = 60;
    var a5 = FG.Alternate(5).AllSubgroups().ToGroupWrapper();
    var allExt = FG.AllExtensions(
            (FG.Abelian(30), FG.Abelian(2)),
            (FG.Abelian(15), FG.Abelian(2, 2)),
            (FG.Abelian(2, 10), FG.Abelian(3))
        )
        .Select(e => e.allSubs.ToGroupWrapper());

    return allExt.AppendIsomorphic([a5], GetAllProds(ord)).Take(GroupExt.A000001[ord]);
}

(AllSubgroups<WElt> subsg, ANameElt[] names)[] AllGroupsFrom33to63()
{
    return 31.Range(33).SelectMany(o => GetAllProds(o))
        .FilterIsomorphic()
        .Concat([..Ord36(), ..Ord40(), ..Ord48(), ..Ord50(), ..Ord54(), ..Ord56(), ..Ord60()])
        .FilterIsomorphic()
        .Naming()
        .ToArray();
}

void AllGroupsUpto63()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();

    var seq = AllGroupsUpTo32()
        .Concat(AllGroupsFrom33to63())
        .CheckMissings()
        .OrderBy(e => e.subsg.Parent.Count())
        .ThenBy(e => e.subsg.Parent.GroupType)
        .ThenBy(e => e.names[0])
        .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
        .ThenBy(e => e.subsg.Infos)
        .ToArray();

    var lvl = Logger.Level;
    Logger.Level = LogLevel.Off;
    foreach (var (subg, names) in seq)
    {
        var g = subg.Parent;
        var rels = Graph.DefiningRelatorsOfGroup(g);
        var relsExpr = rels.Replace(" ", "").Replace(",", ", ").Replace("=", " = ");
        Console.WriteLine($"{g.Count()};{g.Name};{relsExpr}");
    }

    Logger.Level = lvl;
    seq.DisplayNames();

    GlobalStopWatch.Show("End");
    Console.Beep();
}

void Ord64()
{
    Logger.Level = LogLevel.Level1;
    GlobalStopWatch.Restart();

    var g3232 = FG.AllExtensions((FG.Abelian(4, 4), FG.Abelian(2))).Select(e => e.allSubs).Naming()
        .First(e => e.subsg.Infos.ToTuples() == (34, 28, 22)).subsg.Parent.ToCGW();

    var sdps = FG.AllSDPFilterLazy(FG.Abelian(4), FG.Abelian(8))
        .Concat(FG.AllSDPFilterLazy(FG.Abelian(8), FG.Abelian(4)))
        .Select(e => e.AllSubgroups().ToGroupWrapper())
        .FilterIsomorphic()
        .Naming()
        .Select(e => (e.subsg.Parent, FG.Abelian(2)))
        .ToArray();

    var sdp1 = Group.AllSemiDirectProd(FG.Abelian(4, 2), FG.Abelian(4))
        .First(e => e.AllSubgroups().Infos.ToTuples() == (50, 38, 26))
        .ToCGW();
    var sdp2 = Product.Generate(FG.Abelian(2), FG.ModularMaxSdp(4)).ToCGW();
    var listByC2 = sdps.Concat([
        (sdp1, FG.Abelian(2)),
        (sdp2, FG.Abelian(2)),
        (g3232, FG.Abelian(2)),
        (FG.Abelian(8, 4).ToCGW(), FG.Abelian(2)),
        (FG.Abelian(4, 4, 2).ToCGW(), FG.Abelian(2)),
        (FG.Abelian(32).ToCGW(), FG.Abelian(2))
    ]);
    var listByC4_C2C2 = new[]
    {
        (FG.Abelian(8, 2).ToCGW(), FG.Abelian(4)),
        (FG.Abelian(4, 4).ToCGW(), FG.Abelian(4)),
        (FG.Abelian(4, 2, 2).ToCGW(), FG.Abelian(4)),
        (FG.Abelian(8, 2).ToCGW(), FG.Abelian(2, 2))
    };

    var ord8 = new[]
    {
        FG.Abelian(8).ToCGW(),
        FG.Abelian(4, 2).ToCGW(),
        FG.Abelian(2, 2, 2).ToCGW(),
        FG.DihedralSdp(4).ToCGW(),
        FG.Quaternion(8).ToCGW()
    };

    var g1612 = Product.Generate(FG.Quaternion(8), FG.Abelian(2));
    var listSdp64a = FG.MetaCyclicSdp(64).Select(e => e.AllSubgroups().ToGroupWrapper())
        .Concat(ord8.Grid2D(ord8).Select(e => Product.Generate(e.t1, e.t2).AllSubgroups().ToGroupWrapper()))
        .Concat(FG.AllSDPFilter(g1612, FG.Abelian(2, 2)).Select(e => e.AllSubgroups().ToGroupWrapper()))
        .Concat(FG.AllSDPFilter(FG.Abelian(4, 2), FG.DihedralSdp(4)).Select(e => e.AllSubgroups().ToGroupWrapper()))
        .Concat(FG.AllSDPFilter(FG.Abelian(4, 2, 2), FG.Abelian(2, 2)).Select(e => e.AllSubgroups().ToGroupWrapper()));

    var listSdp64b = FG.AllSDPFilter(FG.Abelian(4, 4), FG.Abelian(2)).Select(e => e.AllSubgroups().ToGroupWrapper())
        .Concat(FG.AllSDPFilter(FG.ModularMaxSdp(4), FG.Abelian(2)).Select(e => e.AllSubgroups().ToGroupWrapper()))
        .Append(Product.Generate(FG.Abelian(4), FG.Quaternion(8)).AllSubgroups().ToGroupWrapper())
        .FilterIsomorphic().Naming().ToArray()
        .SelectMany(c => FG.AllSDPFilter(c.subsg.Parent, FG.Abelian(2), trivial: true).Select(e => e.AllSubgroups().ToGroupWrapper()));

    var listAb64 = FG.AllAbelianGroupsOfOrder(64).Select(e => e.AllSubgroups().ToGroupWrapper());
    var listExts = FG.AllExtensions([..listByC2, ..listByC4_C2C2]).Select(e => e.allSubs.ToGroupWrapper());

    var seq = listAb64.Concat([..listSdp64a, ..listSdp64b, ..listExts])
        .FilterIsomorphic()
        .Take(GroupExt.A000001[64])
        .Naming()
        .CheckMissings()
        .OrderBy(e => e.subsg.Parent.Count())
        .ThenBy(e => e.subsg.Parent.GroupType)
        .ThenBy(e => e.names[0])
        .ThenByDescending(e => e.subsg.Parent.ElementsOrders.Values.Max())
        .ThenBy(e => e.subsg.Infos)
        .ToArray();

    var lvl = Logger.Level;
    Logger.Level = LogLevel.Off;
    foreach (var (subg, names) in seq)
    {
        var g = subg.Parent;
        var rels = Graph.DefiningRelatorsOfGroup(g);
        var relsExpr = rels.Replace(" ", "").Replace(",", ", ").Replace("=", " = ");
        Console.WriteLine($"{g.Count()};{g.Name};{relsExpr}");
    }

    Logger.Level = lvl;
    seq.DisplayNames();

    GlobalStopWatch.Show("End");
    Console.Beep();
}

{
    // AllGroupsUpto63();
    // Ord64();
    GroupsFromDB(64);
    // Total Groups:586
    // # End Time:10m28s
}