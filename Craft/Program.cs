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

Dictionary<T, Automorphism<T>> AllInnerAutomorphisms<T>(ConcreteGroup<T> g)
    where T : struct, IElt<T>
{
    var bgAut = new AutomorphismGroup<T>(g);
    var act = Group.ByConjugate(g);
    return g.ToDictionary(
        x => x,
        x => new Automorphism<T>(bgAut,
            Group.IsomorphismMap(g, g, g.GetGenerators().ToDictionary(e => e, e => act(x, e))))
    );
}

(ConcreteGroup<Automorphism<T>> innerAut, Dictionary<T, Automorphism<T>> map) InnerAutomorphismGroup<T>(
    ConcreteGroup<T> g)
    where T : struct, IElt<T>
{
    var allAut = AllInnerAutomorphisms(g);
    var bgAut = allAut.First().Value.AutGroup;
    var autG = Group.Generate($"Inn[{g.Name}]", bgAut, allAut.Values.ToArray());
    return (autG, allAut);
}

Comparer<Perm> CompPerm() => Comparer<Perm>.Create((a, b) =>
{
    var compOrd = a.Order.CompareTo(b.Order);
    if (compOrd != 0)
        return -compOrd;

    var compOrb = a.Orbits.Length.CompareTo(b.Orbits.Length);
    if (compOrb != 0)
        return -compOrb;

    var compFp = a.FixedPoints.Length.CompareTo(b.FixedPoints.Length);
    return compFp;
});

IEnumerable<ConcreteGroup<Perm>> BuildInnerAut(Dictionary<Perm, Automorphism<Perm>> all)
{
    var (a, aut) = all.First();
    var sn = a.Sn;
    var autG = aut.AutGroup;
    var g = autG.Neutral().Domain;
    var inAutG = all.Where(e => e.Value.IsInnerAutomorphism()).ToDictionary();
    var inAutGord = inAutG.Values.Distinct().Count();
    var mapRed = inAutG.Where(e => e.Key.Order != 1 && Group.ElementIsOrder(autG, e.Value, e.Key.Order)).ToDictionary();
    var subGroups = mapRed.Keys.Select(e => new GroupSubset<Perm>([e], Group.GenerateElements(e.Sn, e))).ToHashSet();
    var elts = mapRed.Keys.ToHashSet();
    var allInnerAuts = new HashSet<GroupSubset<Perm>>();
    while (subGroups.Count != 0)
    {
        // Console.WriteLine($"allBags:{allInnerAuts.Count} inAutGord:{inAutGord}");
        // subGroups.Println(s => $"{s.Generators.ToXSet()} {s.Count}", "Before");
        var tmp = subGroups.ToList();
        subGroups.Clear();
        foreach (var set in tmp.OrderBy(s => s.Count))
        {
            var rem = elts.Except(set.Elements).ToHashSet();
            foreach (var e in rem)
            {
                var gens = set.Generators.Append(e).ToHashSet();
                var h = GroupCraft.GenerateElementsLimited(g, set.Elements, gens, inAutGord);
                if (h.Count != 0 && h.All(f => all[f].IsInnerAutomorphism()))
                    subGroups.Add(new(gens, h));
            }
        }

        // subGroups.Println(s => $"{s.Generators.ToXSet()} {s.Count}", "After");
        allInnerAuts.UnionWith(subGroups.Where(b => b.Count == inAutGord));
    }

    // allInnerAuts.Println(b => $"{b.ToXSet()} Is InnerAut:{b.All(e => all[e].IsInnerAutomorphism())}");
    foreach (var bag in allInnerAuts)
    {
        var K = Group.Generate("K", g, bag.ToArray());
        DisplayGroup.HeadElements(K);
        Console.WriteLine();
        yield return K;
    }
}

void InnerAutPerm(ConcreteGroup<Perm> G)
{
    var Z = Group.Zentrum(G);
    DisplayGroup.HeadElements(Z);
    Console.WriteLine($"{Z} IsNormalSubgroup {Group.IsNormalSubgroup(G, Z)}");
    var facts = GroupCraft.AllFactors(G, Z);
    foreach (var (fact, _) in facts)
        DisplayGroup.HeadElements(fact);

    Console.WriteLine();
}

var (success, total, maxTests) = (0, 0, 0);

void SearchAutPerm(ConcreteGroup<Perm> G)
{
    ++total;
    var sn = G.Neutral().Sn;
    var comp = CompPerm();
    var gensG = GroupCraft.RecreateGenerators(G, comp);
    G = Group.Generate(G.Name, G, gensG);
    DisplayGroup.HeadElements(G);
    DisplayGroup.Generators(G);
    var autG = Group.AutomorphismGroup(G);
    DisplayGroup.HeadElements(autG);
    DisplayGroup.Generators(autG);
    Console.WriteLine();
    var act = Group.ByConjugate(sn);
    var setInAut = autG.Where(aut => G.Any(e => aut.AutMap.All(f => act(e, f.Key).Equals(f.Value)))).ToXSet();
    var inAutG = Group.Generate($"Inn[{G.Name}]", autG, setInAut.ToArray());
    
    DisplayGroup.HeadElements(inAutG);
    DisplayGroup.Generators(inAutG);
    var inAutGs = BuildInnerAut(AllInnerAutomorphisms(G)).ToArray();
    // InnerAutPerm(G);
    
    var aut0 = autG.Neutral();
    var gensH = G.Where(e => aut0.AutMap.All(f => act(e, f.Key).Equals(f.Value))).ToArray();
    var Z = Group.Generate($"Z({G.Name})", sn, gensH);
    DisplayGroup.HeadElements(Z);
    var inAutGquo = G.Over(Z);
    
    var cosets = Group.Cosets(G, Z, CosetType.Left);
    cosets.Values.Distinct().Println(l => $"{l} {l.ToXSet()}", $"cosets {G.Name}/{Z.Name}");
    foreach (var K in inAutGs)
        DisplayGroup.AreIsomorphics(K, inAutGquo);
    
    Console.WriteLine();
    Console.WriteLine();
}

void FindAutMetaCyclic(int m, int n, int r)
{
    var (name, sn, gens) = GroupPermutationForm.MetaCyclicGens(m, n, r);
    var G = Group.Generate(name, sn, gens.ToArray());
    GlobalStopWatch.AddLap();
    try
    {
        SearchAutPerm(G);
    }
    catch (Exception e)
    {
        Console.WriteLine(e);
        throw new();
    }

    GlobalStopWatch.Show(G.ShortName);
    Console.WriteLine();
}

void AutPermMetacyclic()
{
    (success, total, maxTests) = (0, 0, 0);
    GlobalStopWatch.AddLap();
    for (int ord = 4; ord <= 32; ord++)
    {
        foreach (var (m, n, r) in GroupPermutationForm.MetaCyclicSdp(ord))
        {
            if ((m, n, r) == (9, 3, 4))
                continue;

            FindAutMetaCyclic(m, n, r);
        }
    }

    GlobalStopWatch.Show($"Missing:{total - success}/{total} Max:{maxTests}");
}

{
    AutPermMetacyclic();
}