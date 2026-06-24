using System.Numerics;
using System.Reflection;
using System.Text;
using Craft;
using Craft.Craft;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.Tools;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

Perm.Style = DisplayPerm.CyclesComma;

string GapExport(Perm[] gens)
{
    var old = Perm.Style;
    Perm.Style = DisplayPerm.Gap;
    var export = $"Group([{gens.Glue(", ")}]);";
    Perm.Style = old;
    return export;
}

(Perm cnD2n, Perm c2D2n, Perm cnAutD2n, Perm[] phiAutD2n, Perm c2nHolD2n) HolD2nPerm(int n)
{
    var pType = IntExt.PrimesDec(n).Select(e => e.Key.Pow(e.Value)).ToArray();
    var a0 = IntExt.PermAndCyclesFromType(pType);
    var n1 = a0.perm.Length;
    var sn = new Sn(2 * n1);
    var bCycles = a0.cycles.SelectMany(e => e.Zip(e.Select(i => i + n1).Reverse().ToArray())).ToArray();
    var cCycles = a0.cycles.Select(o => o.SelectMany(i => new[] { i, i + n1 }).ToArray()).ToArray();
    
    var cnD2n = FG.ConcatPerm(FG.Cycles(n), FG.Cycles(n));
    var c2D2n = sn.OpSeq(bCycles.Select(c => sn.CycleP1([c.First, c.Second])));
    var c2nHolD2n = sn.OpSeq(cCycles.Select(c => sn.CycleP1(c)));

    var cnAutD2n = FG.PaddingRight(FG.Cycles(n), n1);
    var phiAutD2n = FG.AutomorphismDihedralGens(n).gensUn.Select(e => FG.ConcatPerm(e, e)).ToArray();
    
    return (cnD2n, c2D2n, cnAutD2n, phiAutD2n, c2nHolD2n);
}

void RunHolomorphD2n()
{
    GlobalStopWatch.Restart();
    
    for (int n = 3; n <= 32; ++n)
    {
        var (cnD2n, c2D2n, cnAutD2n, phiAutD2n, c2nHolD2n) = HolD2nPerm(n);
        var d2n = Group.Generate($"D{2 * n}", cnD2n.Sn, cnD2n, c2D2n);
        DisplayGroup.HeadOrdersGenerators(d2n);

        var autD2n = Group.Generate($"Aut[{d2n}]", cnD2n.Sn, phiAutD2n.Append(cnAutD2n).ToArray());
        DisplayGroup.HeadOrdersGenerators(autD2n);
        Console.WriteLine($"Aut[D{2 * n}] = C{n} x: {phiAutD2n.Select(e => e.Order).ToAbString().WithParenthesis()}");
        Console.WriteLine();

        var holGens = IntExt.IsPrime(n)
            ? phiAutD2n.Prepend(c2nHolD2n).ToArray() // ord(g1) = 2n and gen of Un ord(g2) = n - 1
            : phiAutD2n.Append(c2D2n).Prepend(c2nHolD2n).ToArray(); // ord(g1) = 2n, ord(g2) = 2 and gens of Un 

        var holD2n = Group.Generate($"Hol[{d2n}]", cnD2n.Sn, holGens);
        DisplayGroup.HeadOrdersGenerators(holD2n);
        Console.WriteLine($"{holD2n} = {d2n} x: {autD2n}");
        Console.WriteLine();

        var d2nNormal = Group.IsNormalSubgroup(holD2n, d2n);
        var autD2nNotNormal = Group.IsNormalSubgroup(holD2n, autD2n);
        var inter = d2n.Intersect(autD2n).Count();
        Console.WriteLine($"{d2n,-10} is normal subgroup of {holD2n,-10} {d2nNormal}");
        Console.WriteLine($"{autD2n,-10} is normal subgroup of {holD2n,-10} {autD2nNotNormal}");
        Console.WriteLine($"{d2n,-10} ∩ {autD2n,10} = {inter}");
        Console.WriteLine();
        Console.WriteLine();
        Console.WriteLine($"hol2:={GapExport(holGens)}");
        Console.WriteLine();
        if (!d2nNormal || autD2nNotNormal || inter != 1)
            throw new();
    }
    
    GlobalStopWatch.Show();
    // GAP
    // d2n:=DihedralGroup(30);autD2n:=AutomorphismGroup(d2n);hol1:=SemidirectProduct(autD2n,d2n);StructureDescription(autD2n);Size(hol1);
    // GeneratorsOfGroup(Image(SmallerDegreePermutationRepresentation(Image(IsomorphismPermGroup(hol1)))));
    // hol2:=Group([(1, 9, 2, 10, 3, 11)(4, 12, 5, 13, 6, 14, 7, 15, 8, 16), (5, 6, 8, 7)(13, 14, 16, 15), (2, 3)(10, 11), (1, 11)(2, 10)(3, 9)(4, 16)(5, 15)(6, 14)(7, 13)(8, 12)]);Size(hol2);
    // IsomorphicSubgroups(hol2,hol1);
}

void RunAutomorphismDihedrals()
{
    GlobalStopWatch.Restart();
    
    for (int n = 3; n <= 32; n++)
    {
        Console.WriteLine($"################################## Aut[D{2 * n}]");
        var (d2n, autD2n) = FG.AutomorphismDihedralPg(n);
        DisplayGroup.HeadOrdersGenerators(d2n);
        DisplayGroup.HeadOrdersGenerators(autD2n); // Aut(D2n) ~ Hol(Cn)
        Console.WriteLine();

        if (!UGCraft.HolCn(FG.Cycles(n)).ToHashSet().SetEquals(autD2n))
            throw new();
    }
    
    GlobalStopWatch.Show();
}

{
    RunAutomorphismDihedrals();
    RunHolomorphD2n();
}