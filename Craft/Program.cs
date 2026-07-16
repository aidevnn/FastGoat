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
    var autbs = new AutomorphismGroup<T>(g);
    var allAut = GroupCraft.AllMorphismsWithPruning(g, g, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(autbs, aut.HomMap))
        .OrderByDescending(aut => Group.GenerateElements(autbs, aut).Count)
        .ToArray();
    var autG = Group.Generate($"Aut[{g.Name}]", autbs, allAut);
    return autG;
}

ConcreteGroup<Perm> RegPermAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var autbs = new AutomorphismGroup<T>(G);
    var sn = new Sn(G.Count());
    var mapEltIdx = G.OrderBy(e => G.ElementsOrders[e]).Index().ToDictionary(e => e.Item, e => e.Index);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(autbs, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return autG;
}

ConcreteGroup<Perm> RegPermZentrumAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var autG = RegPermAutGroup(G);
    var sn = autG.Neutral().Sn;
    var actSn = Group.ByConjugate(sn);
    var zentrumAut = autG.Where(e => autG.GetGenerators().All(aut => actSn.IsInvariant(e, aut))).ToArray();
    var zAutG = Group.Generate($"Z[Aut[{G.Name}]]", sn, zentrumAut);
    var zType = Group.AbelianDecompositions(zAutG).abType.Select(e => e.g).ToArray();
    return Group.Generate($"Z[Aut[{G.Name}]]", sn, zType);
}

ConcreteGroup<Perm> RegPermInnAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var sn = new Sn(G.Count());
    var mapEltIdx = G.OrderBy(e => G.ElementsOrders[e]).Index().ToDictionary(e => e.Item, e => e.Index);
    var act = Group.ByConjugate(G);
    var autbs = Group.AutBase(G);
    var gensInnAut = G.GetGenerators().Select(e => G.GetGenerators().ToDictionary(f => f, f => act(e, f)))
        .Select(pmap => new Automorphism<T>(autbs, Group.IsomorphismMap(G, G, pmap)))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .ToArray();

    return Group.Generate($"InnAut[{G.Name}]", sn, gensInnAut);
}

(ConcreteGroup<Perm> G, ConcreteGroup<Perm> AutG, Dictionary<T, Perm> MapReg)
    RegPermGroupAndAutGroup<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
{
    var sn = new Sn(G.Count());
    var mapEltIdx = G.Index().OrderBy(e => G.ElementsOrders[e.Item])
        .ToDictionary(e => e.Item, e => e.Index);
    var mapIdxElt = mapEltIdx.ToDictionary(e => e.Value, e => e.Key);
    GroupAction<T, Perm> act = (e, p) => sn.CreateElementTable(p.Table.Select(i => mapEltIdx[G.Op(e, mapIdxElt[i])]).ToArray());
    var mapReg = G.ToDictionary(e => e, e => act(e, sn.Neutral()));
    var gensReg = G.GetGenerators().Select(e => mapReg[e]).ToArray();
    var Greg = Group.Generate(G.Name, sn, gensReg);

    var autbs = new AutomorphismGroup<T>(G);
    var allAut = GroupCraft.AllMorphismsWithPruning(G, G, Group.MorphismType.Isomorphism)
        .Select(aut => new Automorphism<T>(autbs, aut.HomMap))
        .Select(e => e.AutMap.ToDictionary(f => mapEltIdx[f.Key], f => mapEltIdx[f.Value]))
        .Select(e => sn.CreateElementTable(e.OrderBy(f => f.Key).Select(f => f.Value).ToArray()))
        .OrderByDescending(e => e.Order)
        .ToArray();
    var autG = Group.Generate($"Aut[{G.Name}]", sn, allAut);
    return (Greg, autG, mapReg);
}

Dictionary<int, (int m, int n, int r)[]> MetaCyclicsM(int m, int maxOrder, bool details = false)
{
    var maxN = maxOrder / m;
    var mCoprimes = IntExt.Coprimes(m).Where(u => u > 1).ToXSet();
    var seqn = (maxN - 1).Range(2).ToDictionary(n => n, n => IntExt.Coprimes(n).ToXSet());
    var all = seqn.SelectMany(e => mCoprimes.Where(u => IntExt.PowMod(u, e.Key, m) == 1)
            .OrderDescending()
            .Select(r => (m, e, r)))
        .DistinctBy(e => (e.e.Key, e.e.Value.Select(i => IntExt.PowMod(e.r, i, m)).ToXSet()))
        .Select(e => (n: e.e.Key, e.r, mCoprimes, nCoprimes: e.e.Value))
        .Where(e => m * e.n <= maxOrder)
        .GroupBy(e => e.r)
        .OrderBy(e => e.Key)
        .ToDictionary(gr => gr.Key, gr => gr.ToArray());

    if (details)
    {
        foreach (var (r, gr) in all)
            gr.Println(e => $"   n:{e.n,-4} {FG.MetaCyclicName(m, e.n, e.r)}", $"m:{m,-4} r:{r}");
    }

    return all.ToDictionary(e => e.Key, e => e.Value.Select(f => (m, f.n, f.r)).ToArray());
}

(string name, int[] abType, int m, int n, int r) RewriteMetaCyclic(int m, int n, int r)
{
    var decM = IntExt.PrimesDec(m).Select(e => e.Key.Pow(e.Value)).ToArray();
    var mAbType = decM.Where(e => r % e == 1).ToArray();
    var m2 = m / mAbType.Aggregate(1, (acc, e) => e * acc);
    var r2 = r % m2;

    var decN = IntExt.PrimesDec(n).Select(e => e.Key.Pow(e.Value)).ToArray();
    var nAbType = decN.Where(e => IntExt.PowMod(r2, e, m2) != 1).ToArray();
    var n2 = n / nAbType.Aggregate(1, (acc, e) => e * acc);
    if (n2 == 1)
        (n2, nAbType) = (n, []);

    var abType = AbelianExt.AbType(mAbType.Concat(nAbType).ToArray()).can;
    var name = FG.MetaCyclicName(m2, n2, r2);
    if (abType[0] != 1)
        name = $"{abType.ToAbString()} x {name}";
    return (name, abType, m2, n2, r2);
}

ConcreteGroup<Perm> HolCn(int n)
{
    var a = FG.Cycles(n);
    return Group.Generate($"Hol[C{n}]", a.Sn, [a, ..FG.UnPermGens(n)]);
}

(string name, ConcreteGroup<Perm> mt, ConcreteGroup<Perm> autMt) AutFrob(int m, int n, int r)
{
    var mt = FG.MetaCyclicPg(m, n, r);
    var (name, abType, m1, n1, r1) = RewriteMetaCyclic(m, n, r);
    var mt1 = FG.MetaCyclicPg(m1, n1, r1);
    var ab = FG.AbelianPerm(abType);
    if (IntExt.Gcd(m1, n1) != 1 || IntExt.Gcd(ab.Order, n1) != 1)
        throw new();

    var autAb = RegPermAutGroup(ab);
    var z1Type = IntExt.Gcd(m, n) == 1 ? Group.AbelianGroupType(RegPermZentrumAutGroup(mt1)) : [1];
    var z1pg = FG.AbelianPerm(z1Type);
    var head = FG.ProductPermGroup(autAb, z1pg);
    if (head.GroupType == GroupType.AbelianGroup)
        head = FG.AbelianPerm(Group.AbelianGroupType(head));

    var tail = HolCn(m1);
    var autMt = head.Order != 1 ? FG.ProductPermGroup([head, tail]) : tail;
    return (name, mt, autMt);
}

void DisplayAutFrob(int m, int n, int r)
{
    var (name, mt, autMt) = AutFrob(m, n, r);
    Console.WriteLine($"##################################### {mt} = {name}");
    DisplayGroup.HeadOrders(mt);
    DisplayGroup.HeadOrders(autMt);
    Console.WriteLine();
}
