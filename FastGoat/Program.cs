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

IEnumerable<AllSubgroups<TableElt>> AllSDP<T1, T2>(ConcreteGroup<T1> N, ConcreteGroup<T2> G)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    Console.WriteLine($"############### AllSDP {N.NameParenthesis()} x: {G.NameParenthesis()}");
    var autG = Group.AutomorphismGroup(G);
    var autN = Group.AutomorphismGroup(N);
    var allOps = Group.AllHomomorphisms(G, autN);
    var ops = allOps.Where(kp => kp.Image().Count() > 1).ToHashSet(new OpByAutEquality<T1, T2>(G, autG, autN));
    Console.WriteLine($"AutG:{autG.Count()} AutN:{autN.Count()}");
    Console.WriteLine($"AllOps:{allOps.Count} Filtered:{ops.Count}");
    foreach (var theta in ops)
        yield return Group.SemiDirectProd(N, theta, G).ToCGTable().AllSubgroups();
}

IEnumerable<AllSubgroups<TableElt>> AllProducts<T1, T2>(ConcreteGroup<T1> G1, ConcreteGroup<T2> G2)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    Console.WriteLine($"G1:{G1.NameParenthesis()} G2:{G2.NameParenthesis()}");
    yield return Product.Generate(G1, G2).ToCGTable().AllSubgroups();

    foreach (var sdp in AllSDP(G1, G2))
        yield return sdp;

    if (G1.Equals(G2))
        yield break;

    foreach (var sdp in AllSDP(G2, G1))
        yield return sdp;
}


string InfosStr(SubGroupsInfos infos) => $"AllSubGr = {infos.AllSubGr}, AllConjsCl = {infos.AllConjsCl}, AllNorms = {infos.AllNorms}";

{
    var fs = File.ReadAllLines("DetailsSmallGroups32.txt");
    var title = "All SmallGroups of order 64";
    var idx = fs.ToList().FindIndex(s => string.Equals(s, title));
    var allInfos = fs.Skip(idx + 1).Chunk(4)
        .GroupBy(e => e[2])
        .ToDictionary(k => k.Key, k => k.Select(e => $"{e[0]}").ToArray());

    var all4 = FG.AllAbelianGroupsOfOrder(4);

    var all32 = FG.AllExtensions(
            (FG.Abelian(16), FG.Abelian(2)),
            (FG.Abelian(2, 8), FG.Abelian(2)),
            (FG.Abelian(4, 4), FG.Abelian(2)),
            (FG.Abelian(2, 4), FG.Abelian(4)),
            (FG.Abelian(2, 4), FG.Abelian(2, 2)))
        .Select(e => e.allSubs.ToTable())
        .Where(e => e.Parent.GroupType == GroupType.NonAbelianGroup || Group.AbelianGroupType(e.Parent).Length < 4)
        .FilterIsomorphic()
        .Naming()
        .Select(e => e.subsg.Parent)
        .ToArray();
    
    var c2 = FG.Abelian(2);
    var ord8 = new[]
    {
        FG.Abelian(8).ToCGTable(),
        FG.Abelian(4, 2).ToCGTable(),
        FG.Abelian(2, 2, 2).ToCGTable(),
        FG.DihedralSdp(4).ToCGTable(),
        FG.Quaternion(8).ToCGTable()
    };
    var ord16 = FG.AllAbelianGroupsOfOrder(16).Where(e => Group.AbelianGroupType(e).Count(c => c == 2) < 4).ToArray();
    
    var list0 = FG.AllAbelianGroupsOfOrder(64).Select(e => e.ToCGTable().AllSubgroups());
    var list1 = ord8.Length.Range().Grid2D().Where(e => e.t1 >= e.t2).Select(e => (t1: ord8[e.t1], t2: ord8[e.t2]))
        .SelectMany(e => AllProducts(e.t1, e.t2));
    var list2 = ord16.Grid2D(all4).SelectMany(e => AllProducts(e.t1, e.t2));
    var list3 = all32.SelectMany(e => AllSDP(e, c2));
    
    var nbOpsMax = 10000;
    var s32 = new Sn(32);
    var g1 = s32[(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16), (17, 18, 19, 20), (21, 22, 23, 24),
        (25, 26, 27, 28), (29, 30, 31, 32)];
    var g2 = s32[(1, 17, 15, 6), (2, 18, 16, 7), (3, 19, 13, 8), (4, 20, 14, 5), (9, 29, 24, 28), (10, 30, 21, 25),
        (11, 31, 22, 26), (12, 32, 23, 27)];
    var g3 = s32[(1, 26, 15, 31), (2, 32, 16, 27), (3, 28, 13, 29), (4, 30, 14, 25), (5, 12, 20, 23),
        (6, 24, 17, 9), (7, 10, 18, 21), (8, 22, 19, 11)];

    var sm3232 = Group.Generate("[(C4 x C4) . C2]pg", s32, g1, g2, g3).ToCGTable();

    var sdps1 = AllProducts(FG.Abelian(4), FG.Abelian(8)).FilterIsomorphic().Naming().Select(e => e.subsg.Parent).ToArray();
    var sdp3 = Product.Generate(c2, FG.ModularMaxSdp(4)).ToCGTable();
    var sdp4 = Group.AllSemiDirectProd(FG.Abelian(2, 4), FG.Abelian(4)).Select(e => e.ToCGTable()).Naming().ToArray()[2];

    var e245 = (nbOpsMax, 0, sm3232, c2); // Group64 Id 172,245
    var tuplesC2 =
        sdps1.Concat([sdp3, sdp4]).Prepend(new Cn(32).ToCGTable()).Select(e => (nbOpsMax, 0, e, c2)).ToArray(); // Group64 Id 11,13,14,22,79,81,82,160,172,180
    var tuplesC4 = new[]
    {
        (nbOpsMax, 0, FG.Abelian(2, 8).ToCGTable(), FG.Abelian(4)),
        (nbOpsMax, 0, FG.Abelian(4, 4).ToCGTable(), FG.Abelian(4))
    }; // Group64 Id 19,22,37,45
    var tuplesC2C2 = new[]
    {
        (nbOpsMax, 0, FG.Abelian(2, 8).ToCGTable(), FG.Abelian(2, 2))
    }; // Group64 Id 49,45,180,160,168,43,172
    
    var list4 = FG.AllExtensions([..tuplesC2, ..tuplesC4, ..tuplesC2C2, e245]).Select(e => e.allSubs.ToTable());

    GlobalStopWatch.Restart();
    var nb = 0;
    
    var all64 = list0.AppendIsomorphic(list1, list2, list3, list4).Take(267).ToArray();
    foreach (var (sub, names) in all64.Naming())
    {
        ++nb;
        sub.Parent.Name = names[0].Name;
        FG.DisplayName(sub.Parent, sub.Infos, names);
        var str = InfosStr(sub.Infos);
        
        var gapNames = allInfos[str];
        if (gapNames.Length == 1)
        {
            Console.WriteLine($"GAP SmallGroup(64) {gapNames[0]}");
            allInfos.Remove(str);
        }
        else
        {
            Console.WriteLine($"Possible GAP SmallGroup(64) {gapNames[0]}");
            if (string.Equals(gapNames[1], "*"))
                allInfos.Remove(str);
            else
                allInfos[str] = gapNames.Skip(1).Append("*").ToArray();
        }

        Console.WriteLine();
    }
    
    Console.WriteLine($"Remaining {allInfos.Sum(e => e.Value.Length)} groups");
    foreach (var (e, f) in allInfos.OrderBy(e => e.Value.Length).ThenBy(e => e.Key))
        f.Println(e);

    Console.WriteLine();
    Console.WriteLine($"Total Groups:{nb}");
    GlobalStopWatch.Show();
    Console.Beep();
    // Remaining 0 groups
    // 
    // Total Groups:267
    // #  Time:2h32m ???????
}
