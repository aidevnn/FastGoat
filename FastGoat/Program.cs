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

void fancyLog()
{
    for (int i = 2; i < 40; i++)
    {
        var cn = new Cn(i);
        foreach (var e in cn.Except([cn.Neutral()]))
            Console.WriteLine("{0} {1} {2}", cn, e, Group.CycleExceptNeutral(cn, e).Glue(" "));
        Console.WriteLine();
    }
}

void Ord16()
{
    GlobalStopWatch.Restart();
    FG.AllExtensions(
            (FG.Abelian(8), FG.Abelian(2)),
            (FG.Abelian(2, 4), FG.Abelian(2)),
            (FG.Abelian(2, 2, 2), FG.Abelian(2)))
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(14)
        .Naming()
        .DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
}

void Ord20_40()
{
    GlobalStopWatch.Restart();
    var all20 = FG.AllExtensions(
            (FG.Abelian(5), FG.Abelian(4)),
            (FG.Abelian(2, 5), FG.Abelian(2)))
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(5)
        .Naming()
        .ToArray();

    var all40 = FG.AllExtensions(
            (FG.Abelian(2, 5), FG.Abelian(4)),
            (FG.Abelian(2, 5), FG.Abelian(2, 2)))
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(14)
        .Naming()
        .ToArray();

    all20.DisplayNames();
    all40.DisplayNames();

    GlobalStopWatch.Show();
    Console.Beep();
}

void Ord32()
{
    GlobalStopWatch.Restart();
    FG.AllExtensions(
            (FG.Abelian(16), FG.Abelian(2)),
            (FG.Abelian(2, 8), FG.Abelian(2)),
            (FG.Abelian(4, 4), FG.Abelian(2)),
            (FG.Abelian(2, 2, 4), FG.Abelian(2)),
            (FG.Abelian(2, 4), FG.Abelian(4)),
            (FG.Abelian(2, 4), FG.Abelian(2, 2)),
            (FG.Abelian(2, 2, 2), FG.Abelian(2, 2)))
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(51)
        .DisplayBox();

    GlobalStopWatch.Show();
    Console.Beep();
    // Total Groups:51
    // #  Time:1m4s
    // 
}

void Ord80()
{
    GlobalStopWatch.Restart();
    FG.AllExtensions(
            (FG.Abelian(2, 5), FG.Abelian(8)),
            (FG.Abelian(2, 2, 5), FG.Abelian(4)),
            (FG.Abelian(4, 5), FG.Abelian(4)),
            (FG.Abelian(4, 5), FG.Abelian(2, 2)),
            (FG.Abelian(2, 2, 5), FG.Abelian(2, 2)),
            (FG.Abelian(2, 2, 2, 2), FG.Abelian(5)))
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(52)
        .DisplayBox();

    GlobalStopWatch.Show("End");
    Console.Beep();
    // Total Groups:52
    // # End Time:6m9s
    // 
}

void Ord81()
{
    GlobalStopWatch.Restart();
    FG.AllExtensions(FG.AllAbelianGroupsOfOrder(27).Select(e => (e, FG.Abelian(3))).ToArray())
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(15)
        .DisplayBox();

    GlobalStopWatch.Show("End");
    Console.Beep();
    // Total Groups:15
    // # End Time:2m26s
}

IEnumerable<AllSubgroups<TableElt>> AllSDP<T1, T2>(ConcreteGroup<T1> N, ConcreteGroup<T2> G)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
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

void ProductOrd32()
{
    GlobalStopWatch.Restart();
    var all2 = FG.AllAbelianGroupsOfOrder(2);
    var all4 = FG.AllAbelianGroupsOfOrder(4);
    var all8 = FG.AllExtensions(
            (FG.Abelian(4), FG.Abelian(2)),
            (FG.Abelian(2, 2), FG.Abelian(2)))
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(5)
        .Naming()
        .DisplayNames()
        .Select(e => e.subsg.Parent)
        .ToArray();

    var all16 = FG.AllExtensions(
            (FG.Abelian(8), FG.Abelian(2)),
            (FG.Abelian(2, 4), FG.Abelian(2)),
            (FG.Abelian(2, 2, 2), FG.Abelian(2)))
        .Select(e => e.allSubs.ToTable())
        .FilterIsomorphic()
        .Take(14)
        .Naming()
        .DisplayNames()
        .Select(e => e.subsg.Parent)
        .ToArray();

    var lt0 = FG.AllAbelianGroupsOfOrder(32).Select(e => e.ToCGTable().AllSubgroups());
    var lt1 = all8.Grid2D(all4).SelectMany(e => AllProducts(e.t1, e.t2));
    var lt2 = all16.Grid2D(all2).SelectMany(e => AllProducts(e.t1, e.t2));
    lt0.AppendIsomorphic(lt1, lt2).DisplayBox();

    GlobalStopWatch.Show();
    Console.Beep();
}
{
    // fancyLog();
    // Ord16();
    // Ord20_40();
    Ord32();
    // Ord80();
    // Ord81();
    TestOrd64();
}
