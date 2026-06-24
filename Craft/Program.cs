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

#region Morphism with pruning

void HomomorphismMap<T1, T2>(IGroup<T1> g, IGroup<T2> h, Dictionary<T1, T2> pMap, Dictionary<T1, T2> map)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    var q = new Queue<T1>();
    q.Enqueue(g.Neutral());
    while (q.Count != 0)
    {
        var e1 = q.Dequeue();
        var f1 = map[e1];
        foreach (var kp in pMap)
        {
            var e2 = kp.Key;
            var f2 = kp.Value;
            var e3 = g.Op(e1, e2);
            var f3 = h.Op(f1, f2);
            if (!map.ContainsKey(e3))
            {
                map[e3] = f3;
                q.Enqueue(e3);
            }
            else if (!map[e3].Equals(f3))
            {
                map.Clear();
                return;
            }
        }
    }
}

IEnumerable<Dictionary<T1, T2>> MorphismPruning<T1, T2>(IGroup<T1> G, IGroup<T2> H, Dictionary<T1, int> sizes,
    Dictionary<T1, T2> pMap, Dictionary<T1, T2> map, (T1 g, T2 a)[][] gens, int idx)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    foreach (var (g, a) in gens[idx])
    {
        var map2 = map.ToDictionary();
        var pMap2 = pMap.ToDictionary();
        pMap2[g] = a;
        HomomorphismMap(G, H, pMap2, map2);
        if (map2.Count == sizes[g])
        {
            if (idx == gens.Length - 1)
                yield return map2;
            else
            {
                foreach (var hom in MorphismPruning(G, H, sizes, pMap2, map2, gens, idx + 1))
                    yield return hom;
            }
        }
    }
}

IEnumerable<Homomorphism<T1, T2>> AllMorphismsWithPruning<T1, T2>(ConcreteGroup<T1> G, ConcreteGroup<T2> H,
    Group.MorphismType mType = Group.MorphismType.Homomorphism)
    where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    var gGens = G.GetGenerators().ToArray();
    bool Filter(int e, int a) => mType == Group.MorphismType.Homomorphism ? e % a == 0 : e == a;
    var hByOrders = H.GroupBy(e => H.ElementsOrders[e])
        .Select(e => (ord: e.Key, elt: e.ToArray()))
        .ToArray();
    var gGensOrders = gGens.Select(e => (g: e, ord: G.ElementsOrders[e])).ToArray();
    var gpMap = gGensOrders.Select(e => hByOrders.Where(a => Filter(e.ord, a.ord)).SelectMany(a => a.elt)
            .Select(a => (e.g, a))
            .ToArray())
        .ToArray();

    var ng = G.Order;
    var sizes = (gGens.Length - 1).SeqLazy(1)
        .ToDictionary(i => gGens[i - 1], i => Group.GenerateElements(G.BaseGroup, gGens.Take(i).ToArray()).Count);
    sizes[gGens.Last()] = ng;
    foreach (var hom in MorphismPruning(G, H, sizes, new(), new() { [G.Neutral()] = H.Neutral() }, gpMap, 0))
    {
        if (mType == Group.MorphismType.Homomorphism)
            yield return new(G, hom);
        else if (hom.Values.Distinct().Count() == ng)
            yield return new(G, hom);
    }
}

bool AreIsomorphic<T1, T2>(ConcreteGroup<T1> G, ConcreteGroup<T2> H) where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    if (G.Order != H.Order || G.GroupType != H.GroupType)
        return false;
    if (!G.ElementsOrdersList().SequenceEqual(H.ElementsOrdersList()))
        return false;

    Homomorphism<T1, T2> nullHom = new();
    var iso = AllMorphismsWithPruning(G, H, Group.MorphismType.Isomorphism).FirstOrDefault(nullHom);
    return !iso.IsNull;
}

#endregion

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