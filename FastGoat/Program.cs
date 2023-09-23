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

// Alejandro Adem, R. James Milgram
// Cohomology of Finite Groups
// I Group Extensions, Simple Algebras and Cohomology
bool CheckTwistedL<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G, Homomorphism<Tg, Automorphism<Tn>> L, MapGroups<Ep2<Tg, Tg>, Tn> w)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    // ω(s, e) = ω(e, s) = e.
    var e = G.Neutral();
    var wes = G.SelectMany(s => new[] { w[new(e, s)], w[new(s, e)] }).ToHashSet();
    if (!wes.SetEquals(new[] { N.Neutral() }))
        return false;

    foreach (var (r, s, t) in G.Grid3D(G, G))
    {
        var lti = L[t].Invert();
        var e0 = N.Op(lti[w[new(r, s)]], w[new(G.Op(r, s), t)]);
        var e1 = N.Op(w[new(s, t)], w[new(r, G.Op(s, t))]);
        if (!e0.Equals(e1))
            return false;
    }

    return true;
}

// Alejandro Adem, R. James Milgram
// Cohomology of Finite Groups
// I Group Extensions, Simple Algebras and Cohomology
HashSet<(Homomorphism<Tg, Automorphism<Tn>>, MapGroups<Ep2<Tg, Tg>, Tn>)> TwistedL<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tg : struct, IElt<Tg> where Tn : struct, IElt<Tn>
{
    var og = G.Count();
    var nG = G.Neutral();
    var nN = N.Neutral();
    var GxG = Product.Generate(G, G);
    var autN = Group.AutomorphismGroup(N);
    var nAN = autN.Neutral();
    var arrG = G.Where(g0 => !g0.Equals(nG)).ToArray();
    var allL = autN.MultiLoop(og - 1).Select(l => arrG.Zip(l).Prepend((nG, nAN)).ToDictionary(e => e.Item1, e => e.Item2))
        .Select(m => new Homomorphism<Tg, Automorphism<Tn>>(G, m)).ToHashSet();
    var allTwisted = new HashSet<(Homomorphism<Tg, Automorphism<Tn>>, MapGroups<Ep2<Tg, Tg>, Tn>)>();
    var we = GxG.Where(e => e.E1.Equals(nG) || e.E2.Equals(nG)).ToArray();
    var we_nN = we.Select(e => (e, nN)).ToArray();
    var rem = GxG.Except(we).ToArray();
    var nb = rem.Length;
    foreach (var l in N.MultiLoop(nb))
    {
        var l0 = l.ToArray();
        var map = new MapGroups<Ep2<Tg, Tg>, Tn>(GxG, N, rem.Zip(l0).Concat(we_nN).ToDictionary(a => a.Item1, a => a.Item2));
        foreach (var L in allL)
        {
            if (CheckTwistedL(N, G, L, map))
                allTwisted.Add((L, map));
        }
    }
    
    return allTwisted;
}

IEnumerable<ConcreteGroup<Ep2<Tn, Tg>>> MapGroups<Tn, Tg>(ConcreteGroup<Tn> N, ConcreteGroup<Tg> G)
    where Tg : struct, IElt<Tg>
    where Tn : struct, IElt<Tn>
{
    var twistedMaps = TwistedL(N, G);
    Console.WriteLine($"{N} \u22b2 E \u2192 {G}; twistedMaps {twistedMaps.Count}");
    var allExts = new List<ConcreteGroup<Ep2<Tn, Tg>>>();
    var count = 0;
    foreach (var (L, map) in twistedMaps)
    {
        var ext = new ExtensionGroup<Tn, Tg>(N, L, map, G);
        var isGroup = Group.IsGroup(ext);
        if (!isGroup)
            continue;

        var extNG = Group.Generate(ext, ext.GetGenerators().ToArray());
        allExts.Add(extNG);

        Console.WriteLine($"{extNG} Count {++count}");
        Console.WriteLine();
        yield return extNG;
    }
}

{
    Console.WriteLine("#############################################");
    var c2 = new Cn(2);
    foreach (var ext in MapGroups(c2, c2))
    {
        DisplayGroup.HeadElements(ext);
        Console.WriteLine($"{ext} is isomorphic to {AbelianInvariantsFactors.Reduce(ext).Glue(" x ", "C{0}")}");
    }
}

{
    var (c2, c4, d8, q8) = (new Cn(2), new Cn(4), FG.Dihedral(4), FG.Quaternion(8));
    Console.WriteLine("#############################################");
    foreach (var ext in MapGroups(c4, c2).ToHashSet(new IsomorphEquality<Ep2<ZnInt, ZnInt>>()))
    {
        DisplayGroup.HeadElements(ext);
        if (ext.GroupType == GroupType.AbelianGroup)
        {
            Console.WriteLine($"{ext} is isomorphic to {AbelianInvariantsFactors.Reduce(ext).Glue(" x ", "C{0}")}");
        }
        else if(ext.ElementsOrders.Values.Order().SequenceEqual(d8.ElementsOrders.Values.Order()))
        {
            Console.WriteLine($"{ext} is isomorphic to D8");
        }
        else if(ext.ElementsOrders.Values.Order().SequenceEqual(q8.ElementsOrders.Values.Order()))
        {
            Console.WriteLine($"{ext} is isomorphic to Q8");
        }
    }
}

{
    Console.WriteLine("#############################################");
    var q8 = FG.Quaternion(8);
    var allExts = MapGroups(FG.DiCyclic(2), new Cn(2)).ToHashSet(new IsomorphEquality<Ep2<Word, ZnInt>>());
    foreach (var ext in allExts)
    {
        var q8ext = Group.IsomorphicsSubgroupsAll(ext, q8);
        if (q8ext.Count != 0)
        {
            Console.WriteLine("#############################################");
            DisplayGroup.HeadElements(ext);
            var e0 = q8ext.First();
            DisplayGroup.HeadElements(e0);
            DisplayGroup.HeadElements(ext.Over(e0));
        }
    }
}

{
    Console.WriteLine("#############################################");
    var s32 = new Sn(32);
    var g1 = s32[(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16), (17, 18, 19, 20), (21, 22, 23, 24),
        (25, 26, 27, 28), (29, 30, 31, 32)];
    var g2 = s32[(1, 17, 15, 6), (2, 18, 16, 7), (3, 19, 13, 8), (4, 20, 14, 5), (9, 29, 24, 28), (10, 30, 21, 25),
        (11, 31, 22, 26), (12, 32, 23, 27)];
    var g3 = s32[(1, 26, 15, 31), (2, 32, 16, 27), (3, 28, 13, 29), (4, 30, 14, 25), (5, 12, 20, 23),
        (6, 24, 17, 9), (7, 10, 18, 21), (8, 22, 19, 11)];

    var sm3232 = Group.Generate("[(C4 x C4) . C2]pg", s32, g1, g2, g3);
    var ext = MapGroups(FG.Abelian(4, 4), new Cn(2)).First(g => g.IsIsomorphicTo(sm3232));
    DisplayGroup.HeadElements(ext);
    DisplayGroup.AreIsomorphics(ext, sm3232);
}

/*
   |(C4 x C4) . C2| = 32
   Type        NonAbelianGroup
   BaseGroup   (C4 x C4) . C2
   
   Elements
   ( 1)[1] = ((0, 0), 0)
   ( 2)[2] = ((0, 2), 0)
   ( 3)[2] = ((2, 0), 0)
   ( 4)[2] = ((2, 2), 0)
   ( 5)[4] = ((0, 0), 1)
   ( 6)[4] = ((0, 1), 0)
   ( 7)[4] = ((0, 1), 1)
   ( 8)[4] = ((0, 2), 1)
   ( 9)[4] = ((0, 3), 0)
   (10)[4] = ((0, 3), 1)
   (11)[4] = ((1, 0), 0)
   (12)[4] = ((1, 0), 1)
   (13)[4] = ((1, 1), 0)
   (14)[4] = ((1, 1), 1)
   (15)[4] = ((1, 2), 0)
   (16)[4] = ((1, 2), 1)
   (17)[4] = ((1, 3), 0)
   (18)[4] = ((1, 3), 1)
   (19)[4] = ((2, 0), 1)
   (20)[4] = ((2, 1), 0)
   (21)[4] = ((2, 1), 1)
   (22)[4] = ((2, 2), 1)
   (23)[4] = ((2, 3), 0)
   (24)[4] = ((2, 3), 1)
   (25)[4] = ((3, 0), 0)
   (26)[4] = ((3, 0), 1)
   (27)[4] = ((3, 1), 0)
   (28)[4] = ((3, 1), 1)
   (29)[4] = ((3, 2), 0)
   (30)[4] = ((3, 2), 1)
   (31)[4] = ((3, 3), 0)
   (32)[4] = ((3, 3), 1)
   
   (C4 x C4) . C2 IsIsomorphicTo [(C4 x C4) . C2]pg : True
*/