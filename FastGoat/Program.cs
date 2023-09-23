using System.Diagnostics;
using System.IO.IsolatedStorage;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
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

IEnumerable<(Homomorphism<T1, T2> i, Homomorphism<T2, T3> p, Homomorphism<T3, T2> s)>
    SplittingGroups<T1, T2, T3>(ConcreteGroup<T1> G, ConcreteGroup<T2> GH, ConcreteGroup<T3> H, bool details = true)
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
    where T3 : struct, IElt<T3>
{
    if (details)
    {
        Console.WriteLine("1 -----> G:{0} -----> GH:{1} -----> H:{2} -----> 1", G, GH, H);
        Console.WriteLine("G  --i--> GH");
        Console.WriteLine("GH --p--> H");
        Console.WriteLine("H  --s--> GH");
        Console.WriteLine("Searching 1 -----> G --i--> GH --p--> H -----> 1 with Im(i)=Ker(p)");
    }

    if (G.Count() * H.Count() != GH.Count())
        yield break;

    var homI = Group.AllHomomorphisms(G, GH);
    var homP = Group.AllHomomorphisms(GH, H);
    var isoS = Group.AllIsomorphisms(H, GH);
    var nullIsos = new Homomorphism<T3, T2>();

    var allImI = homI.Select(hi => (i: hi, im: hi.Image().ToHashSet())).Where(e => e.im.Count > 1).ToArray();
    var allKerP = homP.Select(hp => (p: hp, ker: hp.Kernel())).ToArray();
    var allImIeqKerP = allImI.Grid2D(allKerP)
        .Where(e => e.t1.im.SetEquals(e.t2.ker))
        .Select(e => (e.t1.i, e.t2.p)).ToArray();

    foreach (var ((i, p), k) in allImIeqKerP.Select((e, k) => (e, k)))
    {
        if (details)
        {
            Console.WriteLine($"found {k + 1}");
            Console.WriteLine("Morphism 'i' from G to GH");
            Console.WriteLine("    [{0}]", i);
            Console.WriteLine("Morphism 'p' from GH to H");
            Console.WriteLine("    [{0}]", p);
            Console.WriteLine("Isomorphism 's' from H to GH");
        }

        var found = false;
        foreach (var s in isoS.Where(s0 => s0.Count == H.Count() && H.All(h => p.Domain.Contains(s0[h]) && p[s0[h]].Equals(h))))
        {
            found = true;
            if (details)
                Console.WriteLine("    [{0}]", s);

            yield return (i, p, s);
        }

        if (details)
        {
            if (!found)
            {
                Console.WriteLine("    [not found]");
                Console.WriteLine();
                yield return (i, p, nullIsos);
            }
            else
                Console.WriteLine();
        }
    }
}

// TWISTED ACTIONS AND OBSTRUCTIONS IN GROUP COHOMOLOGY
// IAIN RAEBURN, AIDAN SIMS, AND DANA P. WILLIAMS
bool CheckTwisted<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, MapGroups<Ep2<Tg, Tg>, Tn> w)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    foreach (var (r, s, t) in G.Grid3D(G, G))
    {
        // ω(r, s)ω(rs, t) = ω(s, t)ω(r, st) for r, s, t ∈ G;
        var e0 = N.Op(w[new(r, s)], w[new(G.Op(r, s), t)]);
        var e1 = N.Op(w[new(s, t)], w[new(r, G.Op(s, t))]);
        if (!e0.Equals(e1))
            return false;
    }
    
    // ω(s, e) = ω(e, s) = e.
    var e = G.Neutral();
    var wes = G.SelectMany(s => new[] { w[new(e, s)], w[new(s, e)] }).ToHashSet();
    return wes.Count == 1 && wes.Contains(N.Neutral());
}

HashSet<MapGroups<Ep2<Tg, Tg>, Tn>> Twisted<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tg : struct, IElt<Tg> where Tn : struct, IElt<Tn>
{
    var GxG = Product.Generate(G, G);
    var og2 = GxG.Count();
    var allTwisted = new HashSet<MapGroups<Ep2<Tg, Tg>, Tn>>();
    foreach (var l in N.MultiLoop(og2))
    {
        var l0 = l.ToArray();
        var map = new MapGroups<Ep2<Tg, Tg>, Tn>(GxG, N, GxG.Zip(l0).ToDictionary(a => a.First, a => a.Second));
        if (CheckTwisted(N, G, map))
            allTwisted.Add(map);
    }

    return allTwisted;
}

void MapGroups<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    if (N.GroupType == GroupType.NonAbelianGroup)
        throw new();

    var twistedMaps = Twisted(N, G);
    Console.WriteLine($"{N} \u22b2 E \u2192 {G}; twistedMaps {twistedMaps.Count}");
    var q8 = FG.Quaternion(8);
    var d8 = FG.Dihedral(4);
    foreach (var map in twistedMaps)
    {
        var ext = new ExtensionGroup<Tn, Tg>(N, map, G);
        var isGroup = Group.IsGroup(ext);
        if (!isGroup)
            throw new();
        
        var extNG = Group.Generate(ext, ext.GetGenerators().ToArray());
        DisplayGroup.HeadElements(extNG);
        SplittingGroups(N, extNG, G).First();
        if (extNG.GroupType == GroupType.AbelianGroup)
        {
            Console.WriteLine($"{extNG} is isomorphic to {AbelianInvariantsFactors.Reduce(extNG).Glue(" x ", "C{0}")}");
            Console.WriteLine();
        }
        else if (extNG.Count() == 8 && extNG.IsIsomorphicTo(d8))
        {
            Console.WriteLine($"{extNG} is isomorphic to D8");
            Console.WriteLine();
        }
        else if (extNG.Count() == 8 && extNG.IsIsomorphicTo(q8))
        {
            Console.WriteLine($"{extNG} is isomorphic to Q8");
            Console.WriteLine();
        }

        Console.WriteLine();
    }

    Console.WriteLine();
}

{
    var (c2, c4, v) = (new Cn(2), new Cn(4), FG.Abelian(2, 2));
    MapGroups(c2, c2);
    MapGroups(c4, c2);
    MapGroups(v, c2);
    MapGroups(c2, v);
}