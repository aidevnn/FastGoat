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

ConcreteGroup<Perm> ProductPermGroup(params ConcreteGroup<Perm>[] Gs)
{
    var HnBase = Product.Gp(Gs.Cast<IGroup<Perm>>().ToArray());
    var HnGens = HnBase.GetGenerators().Select(e => FG.ConcatPerm(e.Ei)).ToArray();
    return Group.Generate(HnBase.Name, HnGens[0].Sn, HnGens);
}

bool TypeMatch(ConcreteGroup<Automorphism<Perm>> autG)
{
    var G = autG.Neutral().Domain;
    return autG.GetGenerators().All(aut => G.GetGenerators().All(e => Perm.TypeEquals(e, aut[e])));
}

string GapExport(Perm[] gens)
{
    var old = Perm.Style;
    Perm.Style = DisplayPerm.Gap;
    var export = $"Group([{gens.Glue(", ")}]);";
    Perm.Style = old;
    return export;
}

ConcreteGroup<Automorphism<T>> AutomorphismGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    var bgAut = new AutomorphismGroup<T>(g);
    var allAut = GroupCraft.AllMorphismsWithPruning(g, g, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(bgAut, aut.HomMap))
        .OrderByDescending(aut => Group.GenerateElements(bgAut, aut).Count)
        .ToArray();
    var autG = Group.Generate($"Aut[{g.Name}]", bgAut, allAut);
    return autG;
}

ConcreteGroup<Perm> RegPermAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var bgAut = new AutomorphismGroup<T>(G);
    var sn = new Sn(G.Count());
    var mapEltIdx = G.Index().OrderBy(e => G.ElementsOrders[e.Item]).ToDictionary(e => e.Item, e => e.Index);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(bgAut, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return autG;
}

(ConcreteGroup<Perm> G, ConcreteGroup<Perm> AutG) RegPermGroupAndAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var sn = new Sn(G.Count());
    var mapEltIdx = G.Index().OrderBy(e => G.ElementsOrders[e.Item])
        .ToDictionary(e => e.Item, e => e.Index);
    var mapIdxElt = mapEltIdx.ToDictionary(e => e.Value, e => e.Key);
    GroupAction<T, Perm> act = (e, p) => sn.CreateElementTable(p.Table.Select(i => mapEltIdx[G.Op(e, mapIdxElt[i])]).ToArray());
    var mapReg = G.ToDictionary(e => e, e => act(e, sn.Neutral()));
    var gensReg = G.GetGenerators().Select(e => mapReg[e]).ToArray();
    var Greg = Group.Generate(G.Name, sn, gensReg);
        
    var bgAut = new AutomorphismGroup<T>(G);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(bgAut, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .OrderByDescending(e => e.Order)
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return (Greg, autG);
}

void AutMMnames()
{
    for (int k = 4; k <= 8; k++)
    {
        var mm = FG.ModularMaxPg(k);
        var autMM = RegPermGroupAndAutGroup(mm).AutG;
        DisplayGroup.HeadOrders(autMM);
        var subs = autMM.AllSubgroups();
        subs.DecomposeProducts(subs.ProperNonTrivialNormalSubgroups())
            .Where(e => e.lhs.GroupType == GroupType.AbelianGroup && e.rhs.GroupType == GroupType.AbelianGroup)
            .Select(e => (lhs: Group.AbelianGroupType(e.lhs.Representative),
                rhs: Group.AbelianGroupType(e.rhs.Representative)))
            .DistinctBy(e => e.lhs.ToAbString())
            .Println(e => $"{e.lhs.ToAbString().WithParenthesis()} x: {e.rhs.ToAbString().WithParenthesis()}", autMM.Name);
        Console.WriteLine();
    }
}

void AutQDnames()
{
    for (int k = 4; k <= 7; k++)
    {
        var mm = FG.SemiDihedralPg(k);
        var autMM = RegPermGroupAndAutGroup(mm).AutG;
        DisplayGroup.HeadOrders(autMM);
        var subs = autMM.AllSubgroups();
        subs.DecomposeProducts(subs.ProperNonTrivialNormalSubgroups())
            .Where(e => e.lhs.GroupType == GroupType.AbelianGroup && e.rhs.GroupType == GroupType.AbelianGroup)
            .Select(e => (lhs: Group.AbelianGroupType(e.lhs.Representative),
                rhs: Group.AbelianGroupType(e.rhs.Representative)))
            .DistinctBy(e => e.lhs.ToAbString())
            .Println(e => $"{e.lhs.ToAbString().WithParenthesis()} x: {e.rhs.ToAbString().WithParenthesis()}", autMM.Name);
        Console.WriteLine();
    }
}

void AutMMSDP()
{
    for (int k = 4; k <= 8; k++)
    {
        var n = 2.Pow(k);
        var mm = FG.ModularMaxPg(k);
        var sn = mm.Neutral().Sn;
        DisplayGroup.HeadOrders(mm);
        var autMM = RegPermAutGroup(mm);
        DisplayGroup.HeadOrders(autMM);
        var l = n / 8;
        var (clc2c2, autclc2c2) = RegPermGroupAndAutGroup(FG.AbelianPerm(l, 2, 2));
        var c2c2 = FG.AbelianPerm(2, 2);
        var counter = 0;
        foreach (var theta in GroupCraft.AllMorphismsWithPruning(c2c2, autclc2c2))
        {
            ++counter;
            var gensAutMM = clc2c2.GetGenerators().Concat(c2c2.GetGenerators().Select(e => theta[e])).ToArray();
            var autMMpg = Group.Generate(autMM.Name, sn, gensAutMM);
            if (autMMpg.ElementsOrdersList().SequenceEqual(autMM.ElementsOrdersList()))
            {
                DisplayGroup.HeadOrdersGenerators(autMMpg);
                break;
            }
        }
        
        Console.WriteLine($"try:{counter}");
        Console.WriteLine();
    }
}

void AutQDSDP()
{
    for (int k = 4; k <= 7; k++)
    {
        var n = 2.Pow(k);
        var qd = FG.SemiDihedralPg(k);
        DisplayGroup.HeadOrders(qd);
        var autQD = RegPermAutGroup(qd);
        DisplayGroup.HeadOrders(autQD);
        var l1 = n / 4;
        var l2 = n / 8;
        var (cl, autcl) = RegPermGroupAndAutGroup(FG.AbelianPerm(l1));
        var cl2c2 = FG.AbelianPerm(l2, 2);
        var sn1 = cl.Neutral().Sn;
        var sn2 = cl2c2.Neutral().Sn;
        var counter = 0;
        foreach (var theta in GroupCraft.AllMorphismsWithPruning(cl2c2, autcl))
        {
            ++counter;
            var gensLhs = cl.GetGenerators().Select(e => FG.PaddingRight(e, sn2.N)).ToArray();
            var gensRhs = cl2c2.GetGenerators().Select(e => FG.ConcatPerm(theta[e], e)).ToArray();
            var gensAutQD = gensLhs.Concat(gensRhs).ToArray();
            var autQDpg = Group.Generate(autQD.Name, gensAutQD[0].Sn, gensAutQD);
            if (autQDpg.ElementsOrdersList().SequenceEqual(autQD.ElementsOrdersList()))
            {
                DisplayGroup.HeadOrdersGenerators(autQDpg);
                break;
            }
        }

        Console.WriteLine($"try:{counter}");
        Console.WriteLine();
    }
}

void RunAutMM()
{
    for (int k = 4; k <= 10; k++)
    {
        var n = 2.Pow(k);
        var mm = FG.ModularMaxPg(k);
        DisplayGroup.HeadOrders(mm);
        var l = n / 8;
        var sm = mm.Neutral().Sn;
        var m = sm.N;
        
        var a = sm.OpSeq(4.SeqLazy().Select(i => sm.CycleP1((m / 4).SeqLazy(i, 4).ToArray())));
        var b = sm.OpSeq((m / 2).SeqLazy(0, 2).Select(i => sm.CycleP1([i, i + 1])).ToArray());
        var c = sm.OpSeq(2.SeqLazy().SelectMany(i => (m / 4).SeqLazy(i, 4).Select(j => sm.CycleP1([j, j + 2])))
            .ToArray());
        var d = sm.OpSeq(a.DisjoinCycles.Take(2).Select(e => e ^ (m / 8)));
        
        var sdp = Group.Generate($"(C{l} x C2 x C2) x: C2", sm, a, b, c, d);
        DisplayGroup.HeadOrdersGenerators(sdp);
        var autMM = RegPermAutGroup(mm);
        DisplayGroup.HeadOrders(autMM);
    }
}

void RunAutQD()
{
    for (int k = 4; k <= 8; k++)
    {
        var n = 2.Pow(k);
        var qd = FG.SemiDihedralPg(k);
        DisplayGroup.HeadOrders(qd);
        var l1 = n / 4;
        var l2 = n / 8;
        var (a0, b0, c0) = (FG.Cycles(l1), FG.Cycles(l2), FG.Cycles(2));
        var sn1 = a0.Sn;
        var sn2 = b0.Sn;
        var sn3 = c0.Sn;

        var a = FG.PaddingRight(FG.Cycles(l1), sn2.N + sn3.N);
        var b = FG.Padding(sn1.N, b0, sn3.N);
        var c1 = sn1.OpSeq((l1 / 2 - 1).SeqLazy(1).Select(i => sn1.CycleP1([i, l1 - i])));
        var c = FG.ConcatPerm(FG.PaddingRight(c1, sn2.N), c0);
        var sdp = Group.Generate($"C{l1} x: (C{l2} x C2)", a.Sn, a, b, c);
        DisplayGroup.HeadOrders(sdp);
        var autQD = RegPermAutGroup(qd);
        DisplayGroup.HeadOrders(autQD);
    }
}

void TestMM()
{
    AutMMnames();
    AutMMSDP();
    RunAutMM();
}

void TestQD()
{
    AutQDnames();
    AutQDSDP();
    RunAutQD();
}