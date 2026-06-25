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

(ConcreteGroup<Perm> Dicm, ConcreteGroup<Perm> AutDicm) AutomorphismDicyclic(int m)
{
    if (int.IsPow2(m))
    {
        var Qm = FG.QuaternionPg(m * 4);
        var autQm = m == 2 ? FG.Symmetric(4) : FG.AutomorphismDihedralPg(m * 2).AutD2n;
        autQm.Name = $"Aut[{Qm}]";
        return (Qm, autQm);
    }
    else if (m % 2 == 1)
    {
        var Dicm = FG.DicyclicPg(m);
        var (a0, b0, gensUn0) = FG.AutomorphismDihedralGens(m);
        var a = FG.PaddingRight(a0, 2);
        var b = FG.PaddingRight(b0, 2);
        var c = FG.PaddingLeft(FG.Cycles(2), b0.Sn.N);
        var gensUn = gensUn0.Select(e => FG.PaddingRight(e, 2)).ToArray();
        var AutDicm = Group.Generate($"Aut[{Dicm}]", a.Sn, gensUn.Concat([a, b, c]).ToArray());
        return (Dicm, AutDicm);
    }
    else
    {
        var Dicm = FG.DicyclicPg(m);
        var k = IntExt.PrimesDecomposition(m).Count(i => i == 2);
        var n = m / (1 << k);

        var (a0, b0, gensUn0) = FG.AutomorphismDihedralGens(n);
        var autD2n = Group.Generate($"Aut[D{2 * n}]", a0.Sn, gensUn0.Prepend(a0).ToArray());

        var q = 1 << (k + 1);
        var (a1, b1, gensUn1) = FG.AutomorphismDihedralGens(q);
        var autD2q = Group.Generate($"Aut[D{2 * q}]", a1.Sn, gensUn1.Prepend(a1).ToArray());
    
        var hom = gensUn1.ToDictionary(e => e, _ => a0.Sn.Neutral());
        hom[a1] = Group.GenerateElements(a0.Sn, gensUn0).First(e => e.Order == 2);

        var gens1 = autD2n.GetGenerators().Select(c => FG.PaddingRight(c, a1.Sn.N)).ToArray();
        var gens2 = autD2q.GetGenerators().Select(c => FG.ConcatPerm(hom[c], c)).ToArray();
        var sn = gens1[0].Sn;
        var gensAutDicm = gens1.Concat(gens2).ToArray();
        var AutDicm = Group.Generate($"Aut[{Dicm}]", sn, gensAutDicm);
        return (Dicm, AutDicm);
    }
}

{
    for (int m = 2; m <= 64; m++)
    {
        var (dicm, autDicm) = AutomorphismDicyclic(m);
        DisplayGroup.HeadOrders(dicm);
        DisplayGroup.HeadOrdersGenerators(autDicm);
        if (m <= 16)
        {
            var autDicmReg = FG.RegPermAutGroup(dicm).autG;
            if (!GroupCraft.AreIsomorphic(autDicmReg, autDicm))
                throw new();
        }

        Console.WriteLine($"autDic{m}:={GapExport(autDicm.GetGenerators().ToArray())}");
        Console.WriteLine();
    }
}
