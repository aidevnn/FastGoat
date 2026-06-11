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
    Group.MorphismType mType)
    where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    var gGens = G.GetGenerators().ToArray();
    bool Filter(int e, int a) => mType == Group.MorphismType.Homomorphism ? e % a == 0 : e == a;
    var g2ByOrders = H.GroupBy(e => H.ElementsOrders[e])
        .Select(e => (ord: e.Key, elt: e.ToArray()))
        .ToArray();
    var gGensOrders = gGens.Select(e => (g: e, ord: G.ElementsOrders[e])).ToArray();
    var gpMap = gGensOrders.Select(e => g2ByOrders.Where(a => Filter(e.ord, a.ord)).SelectMany(a => a.elt)
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

bool AutomorphismReconstruction(ConcreteGroup<Perm> G)
{
    var autG = Group.AutomorphismGroup(G);
    var autGgens = GroupCraft.RecreateGenerators(autG);
    if (autGgens.Any(aut => G.GetGenerators().Any(e => !Perm.TypeEquals(e, aut[e]))))
    {
        // Console.WriteLine($"{autG.Name} not found in {G.Neutral().Sn} TypeMisMatch");
        return false;
    }

    DisplayGroup.HeadOrdersGenerators(G);
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
    Dictionary<Automorphism<Perm>, Perm> emptyMap = new();
    Dictionary<Automorphism<Perm>, Perm> pmap = new();
    Dictionary<Automorphism<Perm>, Perm> map = new() { [autG.Neutral()] = G.Neutral() };
    var iso = MorphismPruning(autG, G.BaseGroup, sizes, pmap, map, gensTuples, 0).FirstOrDefault(emptyMap);
    if (iso.Count == 0)
    {
        // Console.WriteLine($"{autG.Name} not found in {G.Neutral().Sn} Isomorphism");
        return false;
    }

    var autGperm = Group.Generate(autG.Name + "pg", G.BaseGroup, autGgens.Select(aut => iso[aut]).ToArray());
    DisplayGroup.HeadOrdersGenerators(autGperm);
    Console.WriteLine($"{autG} Is Isomorphic To {autGperm}");
    Console.WriteLine();
    return true;
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

{
    GlobalStopWatch.Restart();
    for (int m = 2; m <= 16; m++)
    {
        if (!int.IsPow2(m))
            AutomorphismReconstruction(FG.DicyclicPg(m));
    }
    GlobalStopWatch.Show();
}