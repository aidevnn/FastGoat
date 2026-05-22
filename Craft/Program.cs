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

IEnumerable<Homomorphism<T1, T2>> HomomorphismFromKernel<T1, T2>(ConcreteGroup<T1> g1, ConcreteGroup<T2> g2, HashSet<T1> kerGens)
    where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    if (!kerGens.IsSubsetOf(g1))
        throw new();
    
    var kerMap = kerGens.Select(e => (g: e, a: g2.Neutral())).ToArray();
    var Ker = Group.GenerateElements(g1, kerGens.ToArray());
    var gGens = g1.GetGenerators().Except(Ker).ToArray();
    
    var g2ByOrders = g2.GroupBy(e => g2.ElementsOrders[e])
        .Select(e => (ord: e.Key, elt: e.ToArray()))
        .ToArray();
    
    var gGensOrders = gGens.Select(e => (g: e, ord: g1.ElementsOrders[e])).ToArray();
    var gpMap = gGensOrders.Select(e => g2ByOrders.Where(a => e.ord % a.ord == 0).SelectMany(a => a.elt)
            .Select(a => (e.g, a)).ToArray()).ToArray();

    var ng = g1.Count();
    foreach (var arr in gpMap.MultiLoop())
    {
        var map = arr.Concat(kerMap).ToDictionary(t => t.g, t => t.a);
        var hom = Group.HomomorphismMap(g1, g2, map);
        if (hom.Count == ng)
            yield return new(g1, hom);
    }
}

Dictionary<ConcreteGroup<T>, XSet<XSet<T>>> NormalsSubgroupsAndFactors<T>(ConcreteGroup<T> G, bool proper = false)
    where T : struct, IElt<T>
{
    return Group.AllHomomorphisms(G, G).GroupBy(hom => hom.Kernel().ToXSet())
        .ToDictionary(e => e.Key, e => e.Select(hom => hom.Image().ToXSet()).ToXSet())
        .Where(e => !proper || (e.Key.Count > 1 && e.Key.Count < G.Order))
        .OrderByDescending(e => e.Key)
        .ToDictionary(e => Group.Generate("N", G, e.Key.ToArray()), e => e.Value);
}

IEnumerable<ConcreteGroup<T>> Factor<T>(ConcreteGroup<T> G, ConcreteGroup<T> H) where T : struct, IElt<T>
{
    var ordK = G.Order / H.Order;
    return HomomorphismFromKernel(G, G, H.GetGenerators().ToHashSet())
        .Where(hom => !hom.IsNull)
        .Select(hom => hom.Image().ToXSet())
        .Where(img => img.Count == ordK && img.X.Intersect(H).Count() == 1)
        .Distinct()
        .Select(img => Group.Generate("K", G, img.ToArray()));
}

{
    var G = FG.SemiDihedral(4);
    DisplayGroup.HeadElements(G);
    foreach (var (N, Factors) in NormalsSubgroupsAndFactors(G, proper: true))
    {
        Console.WriteLine("############################");
        DisplayGroup.HeadElements(N);
        foreach (var K in Factor(G, N).Take(1))
        {
            DisplayGroup.HeadElements(K);
            Console.WriteLine("K is Normal in G:{0}", Group.IsNormalSubgroup(G, K));
            Console.WriteLine("K is in Factors:{0}", Factors.Contains(K.ToXSet()));
            Console.WriteLine();
        }

        Console.WriteLine($"Factors:{Factors.Count} Disjoint:{Factors.Count(f => f.Intersect(N).Count() == 1)}");
        Console.WriteLine();
    }
}