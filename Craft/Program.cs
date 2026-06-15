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

void TestMorphPruning()
{
    for (int ord = 1; ord <= 32; ord++)
    {
        foreach (var G0 in FG.AllGroupsOfOrder(ord).Where(g => g.GroupType == GroupType.NonAbelianGroup))
        {
            var G = G0.ToPermGroup().Item1;
            Console.WriteLine(G.ShortName);
            if (ord < 32)
            {
                var mType = Group.MorphismType.Homomorphism;
                var allHoms1 = AllMorphismsWithPruning(G, G, mType).ToArray();
                Console.WriteLine($"allHoms1:{allHoms1.Length}");
                var allHoms2 = Group.AllMorphisms(G, G, mType).ToArray();
                Console.WriteLine($"allHoms2:{allHoms2.Length}");
                var test = allHoms1.Order().SequenceEqual(allHoms2.Order());
                Console.WriteLine($"seq equal:{test}");
                if (!test)
                    throw new();
                Console.WriteLine();
            }

            {
                var mType = Group.MorphismType.Isomorphism;
                var allIsos1 = AllMorphismsWithPruning(G, G, mType).ToArray();
                Console.WriteLine($"allIsos1:{allIsos1.Length}");
                var allIsos2 = Group.AllMorphisms(G, G, mType).ToArray();
                Console.WriteLine($"allIsos2:{allIsos2.Length}");
                var test = allIsos1.Order().SequenceEqual(allIsos2.Order());
                Console.WriteLine($"seq equal:{test}");
                Console.WriteLine();
                if (!test)
                    throw new();
                Console.WriteLine();
            }
        }
    }
}

void BenchMorphPruning()
{
    var mType = Group.MorphismType.Isomorphism;
    foreach (var G0 in FG.AllGroupsOfOrder(16).Where(g => g.GroupType == GroupType.NonAbelianGroup))
    {
        var (G, autG, _) = FG.RegPermAutGroup(G0);
        if (autG.GetGenerators().Count() != 4)
            continue;
        Console.WriteLine(G.ShortName);
        Console.WriteLine(autG.ShortName);

        GlobalStopWatch.Bench(5, "Old All", () => Group.AllMorphisms(autG, autG, mType).ToArray());
        GlobalStopWatch.Bench(5, "New All", () => AllMorphismsWithPruning(autG, autG, mType).ToArray());
        Console.WriteLine();
    }
}

(Homomorphism<T, Automorphism<T>> opsByAut, ConcreteGroup<Automorphism<T>> innAutG)
    InnerAutomorphismGroup<T>(ConcreteGroup<T> g, ConcreteGroup<T> autG) where T : struct, IElt<T>
{
    var bgAut = new AutomorphismGroup<T>(g);
    var act = Group.ByConjugate(autG);
    var allAut = autG.ToDictionary(x => x,
        x => new Automorphism<T>(bgAut,
            Group.IsomorphismMap(g, g, g.GetGenerators().ToDictionary(e => e, e => act(x, e)))));
    var innAutG = Group.Generate($"Inn[{autG.Name}]", bgAut, allAut.Values.ToArray());
    return (new Homomorphism<T, Automorphism<T>>(autG, allAut), innAutG);
}

(ConcreteGroup<Perm> G, ConcreteGroup<Perm> AutG)
    AutomorphismPerm(ConcreteGroup<Perm> G, ConcreteGroup<Automorphism<Perm>> autG)
{
    var sn = G.Neutral().Sn;
    if (!TypeMatch(autG))
        return (sn.C1, sn.C1);

    if (G.Order == sn.N)
    {
        Console.WriteLine($"Regular permutation {G}");
        var (G0, autG0, _) = FG.RegPermAutGroup(G);
        return (G0, autG0);
    }

    var autGgens = GroupCraft.RecreateGenerators(autG);
    var autGensImages =
        autGgens.ToDictionary(aut => aut,
            aut => UGCraft.InnerAutMatrix(aut).Where(e => e.Order == autG.ElementsOrders[aut]).Select(a => (g: aut, a))
                .ToArray());
    var ng = autG.Order;
    var sizes = (autGgens.Length - 1).SeqLazy(1)
        .ToDictionary(i => autGgens[i - 1],
            i => Group.GenerateElements(autG.BaseGroup, autGgens.Take(i).ToArray()).Count);
    sizes[autGgens.Last()] = ng;
    var gensTuples = autGgens.Select(aut => autGensImages[aut]).ToArray();
    Dictionary<Automorphism<Perm>, Perm> pmap = new();
    Dictionary<Automorphism<Perm>, Perm> map = new() { [autG.Neutral()] = G.Neutral() };
    foreach (var iso in MorphismPruning(autG, G.BaseGroup, sizes, pmap, map, gensTuples, 0))
    {
        var autGgensPg = autGgens.Select(aut => iso[aut]).ToArray();
        var autGperm = Group.Generate(autG.Name + "pg", G.BaseGroup, autGgensPg);
        if (AreIsomorphic(autG, autGperm))
        {
            Console.WriteLine($"{autG} Is Isomorphic To {autGperm} in {G.BaseGroup}");
            return (G, autGperm);
        }
    }

    return (G.C1, G.C1);
}

bool TypeMatch(ConcreteGroup<Automorphism<Perm>> autG)
{
    var G = autG.Neutral().Domain;
    return autG.GetGenerators().All(aut => G.GetGenerators().All(e => Perm.TypeEquals(e, aut[e])));
}

void RunTests()
{
    TestMorphPruning();
    // |(C3 x C3) x: C3| = 27
    // allHoms1:729
    // allHoms2:729
    // seq equal:True
    // 
    // allIsos1:432
    // allIsos2:432
    // seq equal:True
    // 

    BenchMorphPruning();
    // |Q16| = 16
    // |Aut[Q16]| = 32
    // # Old All Avg Time:388 ms Dev:11.923
    // # New All Avg Time:7 ms Dev:3.382
    // 
    // |M(4x:4)3| = 16
    // |Aut[M(4x:4)3]| = 32
    // # Old All Avg Time:1451 ms Dev:38.436
    // # New All Avg Time:483 ms Dev:7.574
    // 
    // |C2 x D8| = 16
    // |Aut[C2 x D8]| = 64
    // # Old All Avg Time:15817 ms Dev:1579.559
    // # New All Avg Time:234 ms Dev:12.986
    // 
    // |(C2 x C2) x: C4| = 16
    // |Aut[(C2 x C2) x: C4]| = 32
    // # Old All Avg Time:3192 ms Dev:67.689
    // # New All Avg Time:1035 ms Dev:21.479
    // 
}

IEnumerable<(Perm a, Perm b, int abCount)> rUmToPerm(int m, int n, int r, int dim = 1)
{
    var a = FG.ConcatPerm(Enumerable.Repeat(FG.Cycles(m), dim).ToArray());
    var orbx = a.Orbits.Select(o => o.ToXSet()).ToArray();
    var idxOrbx = orbx.SelectMany(o => o.Select(i => (i, o))).ToDictionary(e => e.i, e => e.o);
    foreach (var b in UGCraft.InnerAut(a, a ^ r).Where(b => b.Order != 1 && n % b.Order == 0))
    {
        var ab = GroupCraft.GenerateElementsLimited(a.Sn, [a, b], m * n);
        if (m * n % ab.Count == 0 && b.Orbits.All(o => o.Select(i => idxOrbx[i]).Distinct().Count() == dim))
            yield return (a, b, ab.Count);
    }
}

IEnumerable<(Perm a, Perm b, ConcreteGroup<Perm> mt)> MetaCyclicPg(int m, int n, int r, int dim = 1)
{
    var mt1 = FG.MetaCyclicPg(m, n, r);
    foreach (var (a0, b0, count) in rUmToPerm(m, n, r, dim).Order().DistinctBy(e => e.b.PermTypeStr))
    {
        if (count == m * n)
        {
            var mt2 = Group.Generate(mt1.Name, a0.Sn, a0, b0);
            if (AreIsomorphic(mt1, mt2))
                yield return (a0, b0, mt2);
            else
                Console.WriteLine("mtGensProblems1");
        }
        else
        {
            var n0 = IntExt.DividorsInt(n).Order().First(x => IntExt.Lcm(x, b0.Order) == n);
            var b1 = FG.Cycles(n0);
            var a = FG.PaddingRight(a0, b1.Sn.N);
            var b = FG.ConcatPerm(b0, b1);
            var mt2 = Group.Generate(mt1.Name, a.Sn, a, b);
            if (AreIsomorphic(mt1, mt2))
                yield return (a, b, mt2);
            else if (n0 != n && mt2.Order < m * n)
            {
                var b2 = FG.Cycles(n);
                var a3 = FG.PaddingRight(a0, b2.Sn.N);
                var b3 = FG.ConcatPerm(b0, b2);
                var mt3 = Group.Generate(mt1.Name, a3.Sn, a3, b3);
                if (AreIsomorphic(mt1, mt3))
                    yield return (a3, b3, mt3);
                else
                    Console.WriteLine("mtGensProblems2");
            }
        }
    }
}

Perm[] InnerHol(ConcreteGroup<Perm> G, ConcreteGroup<Perm> autG)
{
    var holOrd = G.Order * autG.Order;
    var holGens = autG.GetGenerators().Concat(G.GetGenerators()).ToArray();
    var sn = holGens[0].Sn;
    var hol = GroupCraft.GenerateElementsLimited(sn, holGens, holOrd);
    if (hol.Count == holOrd)
        return holGens;
    else
        return [];
}

Perm[] OuterHol(ConcreteGroup<Perm> G, ConcreteGroup<Perm> autG)
{
    var holOrd = G.Order * autG.Order;
    var (theta, _) = InnerAutomorphismGroup(G, autG);
    var gensConcats = autG.GetGenerators()
        .SelectMany(g => G.GetGenerators().Select(a => FG.ConcatPerm(theta[g][a], g)))
        .ToArray();
    var sn = gensConcats[0].Sn;
    var gensPads = G.GetGenerators().Select(g => FG.PaddingRight(g, sn.N - g.Sn.N))
        .ToArray();
    var holGens = gensPads.Concat(gensConcats).ToArray();
    var hol = GroupCraft.GenerateElementsLimited(sn, holGens, holOrd);
    if (hol.Count == holOrd)
        return holGens;
    else
        return [];
}

// void RunHolomorphsMetaCyclic()
{
    var stats = new List<string>();
    GlobalStopWatch.Restart();
    for (int ord = 1; ord <= 64; ord++)
    {
        foreach (var (m, n, r) in FG.FrobeniusCoefs(ord))
        {
            var name = FG.MetaCyclicName(m, n, r);
            Console.WriteLine($"{name}");
            if (IntExt.PrimesDec(m * n).Count > 1)
            {
                var dims = new Dictionary<(string, Sn), Perm[]>();
                foreach (var dim in new[] { 1, 2 })
                {
                    foreach (var (a, b, M1) in MetaCyclicPg(m, n, r, dim))
                    {
                        var autM1 = Group.AutomorphismGroup(M1);
                        var (M2, autM2pg) = AutomorphismPerm(M1, autM1);
                        if (autM2pg.Order == autM1.Order)
                        {
                            if (dim == 1)
                            {
                                var holGens = OuterHol(M2, autM2pg); 
                                if(holGens.Length > 1)
                                    dims[("outer", holGens[0].Sn)] = holGens;
                            }
                            
                            if (dim == 2)
                            {
                                var holGens = InnerHol(M2, autM2pg);
                                if(holGens.Length > 1)
                                    dims[("inner", holGens[0].Sn)] = holGens;
                            }
                        }
                    }
                }

                var ((case0, snBest), holG) = dims.MinBy(e => e.Key.Item2.N);
                holG.Println($"Hol[{name}] in {snBest} {case0}");
                stats.Add(case0);
                dims.Keys.Println("Dims");
            }
            else
            {
                var (M2, autM2, _) = FG.RegPermAutGroup(FG.MetaCyclicPg(m, n, r));
                var holOrd = M2.Order * autM2.Order;
                var holGens = autM2.GetGenerators().Concat(M2.GetGenerators()).ToArray();
                var sn = holGens[0].Sn;
                var hol = GroupCraft.GenerateElementsLimited(sn, holGens, holOrd);
                if (hol.Count == holOrd)
                {
                    holGens.Println($"Hol[{M2}] in {sn}");
                    stats.Add("inner");
                }
                else
                    Console.WriteLine($"    Problem Hol[{M2}] Reg");
            }
        }
    }

    GlobalStopWatch.Show();
    stats.GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count()).Println("Stats");
}
