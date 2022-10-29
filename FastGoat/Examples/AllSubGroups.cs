using FastGoat.Structures.CartesianProduct;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace FastGoat.Examples;

public static class AllSubGroups
{
    static (HashSet<HashSet<T>> im, HashSet<HashSet<T>> ker) AllImKer<T>(ConcreteGroup<T> g0, ConcreteGroup<T> g1)
        where T : struct, IElt<T>
    {
        HashSet<HashSet<T>> all = new HashSet<HashSet<T>>(new SetEquality<T>());
        HashSet<HashSet<T>> allIm = new HashSet<HashSet<T>>(new SetEquality<T>());
        HashSet<HashSet<T>> allKer = new HashSet<HashSet<T>>(new SetEquality<T>());
        var hom = Group.AllHomomorphisms(g0, g1);
        foreach (var m in hom)
        {
            var im = m.Image().ToHashSet();
            var ker = m.Kernel().ToHashSet();
            allIm.Add(im);
            allKer.Add(ker);
        }

        return (allIm, allKer);
    }

    static void AllNormalSubGroups<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var allHoms = Group.AllHomomorphisms(g, g);
        var allKer = allHoms.Select(p => p.Kernel().ToHashSet())
            .ToHashSet(new SetEquality<T>());
        foreach (var ker in allKer.OrderBy(a0 => a0.Count))
        {
            var gKer = Group.Generate($"Ker{ker.Count}", g, ker.ToArray());
            DisplayGroup.Head(g.Over(gKer));
        }
    }

    public static void PrimeCycleSimple()
    {
        var g = new Cn(11);
        AllNormalSubGroups(g);
    }

    public static void PrimeCycleSimpleAnother()
    {
        var g = new Cn(53);
        AllNormalSubGroups(g);
    }

    public static void CyclicNotSimple()
    {
        var g = new Cn(51);
        AllNormalSubGroups(g);
    }

    public static void CyclicNotSimpleAnother()
    {
        var g = Product.Generate(new Cn(3), new Cn(4));
        AllNormalSubGroups(g);
    }

    public static void A4NotSimple()
    {
        var s4 = new Sn(4);
        var g = Group.Generate("A4", s4, s4[(1, 3), (2, 4)], s4[(1, 2, 3)]);
        AllNormalSubGroups(g);
    }

    public static void A5Simple()
    {
        var s5 = new Sn(5);
        var g = Group.Generate("A5", s5, s5[(1, 3, 5)], s5[(1, 2, 3, 4, 5)]);
        AllNormalSubGroups(g);
    }

    public static void Symm6NormalSubGroups()
    {
        var g = Group.Create(new Sn(6));
        AllNormalSubGroups(g);
    }

    public static void AbelianSubGroups()
    {
        var g = Product.Generate(new Cn(4), new Cn(8));
        var (allIm, allKer) = AllImKer(g, g);
        var all = allIm.Union(allKer).ToHashSet(new SetEquality<Ep2<ZnInt, ZnInt>>());
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void AbelianSubGroupsAnother()
    {
        var g = Product.Generate(new Cn(6), new Cn(9));
        var (allIm, allKer) = AllImKer(g, g);
        var all = allIm.Union(allKer).ToHashSet(new SetEquality<Ep2<ZnInt, ZnInt>>());
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void AbelianSubGroupsEncore()
    {
        var g = Product.Generate(new Cn(2), new Cn(4), new Cn(8));
        var (allIm, allKer) = AllImKer(g, g);
        var all = allIm.Union(allKer).ToHashSet(new SetEquality<Ep3<ZnInt, ZnInt, ZnInt>>());
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void Symm3SubGroups()
    {
        var g = Group.Create(new Sn(3));
        var all = new HashSet<HashSet<Perm>>(new SetEquality<Perm>());
        all.UnionWith(g.Select(e => Group.Generate(g, e).ToHashSet()));
        Console.WriteLine($"Monogenics SubGroups from elements : {all.Count()}");

        var (_, allKer) = AllImKer(g, g);
        all.UnionWith(allKer);
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void Symm4SubGroups()
    {
        var g = Group.Create(new Sn(4));
        var all = new HashSet<HashSet<Perm>>(new SetEquality<Perm>());
        all.UnionWith(g.Select(e => Group.Generate(g, e).ToHashSet()));
        Console.WriteLine($"Monogenics SubGroups from elements : {all.Count()}");

        var (allIm, allKer) = AllImKer(g, g);
        all.UnionWith(allIm);
        all.UnionWith(allKer);
        Console.WriteLine($"All Monogenics and Normals SubGroups : {all.Count()}");

        // Dihedral is added by hand
        var d8 = Group.Generate(g, g[(1, 4), (2, 3)], g[(1, 2, 3, 4)]);
        var (allIm8, allKer8) = AllImKer(d8, g);
        all.UnionWith(allIm8);
        all.UnionWith(allKer8);
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void Symm5SubGroups()
    {
        var g = Group.Create(new Sn(5));
        var all = new HashSet<HashSet<Perm>>(new SetEquality<Perm>());
        all.UnionWith(g.Select(e => Group.Generate(g, e).ToHashSet()));
        Console.WriteLine($"Monogenics SubGroups from elements : {all.Count()}");

        var (allIm, allKer) = AllImKer(g, g);
        all.UnionWith(allIm);
        all.UnionWith(allKer);
        Console.WriteLine($"All Monogenics and Normals SubGroups : {all.Count()}");

        // S3 subgroups and isomorphic are added by hand
        var s3 = Group.Generate(g, g[(1, 2)], g[(1, 2, 3)]);
        var (allIm6, allKer6) = AllImKer(s3, g);
        all.UnionWith(allIm6);
        all.UnionWith(allKer6);
        Console.WriteLine($"With S3 subgroups : {all.Count()}");

        // Dihedral 8 subgroups and isomorphic are by hand
        var d8 = Group.Generate(g, g[(1, 4), (2, 3)], g[(1, 2, 3, 4)]);
        var (allIm8, allKer8) = AllImKer(d8, g);
        all.UnionWith(allIm8);
        all.UnionWith(allKer8);
        Console.WriteLine($"With D8 subgroups : {all.Count()}");

        // Dihedral 10 subgroups and isomorphic are by hand
        var d10 = Group.Generate(g, g[(2, 5), (3, 4)], g[(1, 2, 3, 4, 5)]);
        var (allImD10, allKerD10) = AllImKer(d10, g);
        all.UnionWith(allImD10);
        all.UnionWith(allKerD10);
        Console.WriteLine($"With D10 subgroups : {all.Count()}");

        // A4 subgroups and isomorphic are added by hand
        var a4 = Group.Generate("A4", g, g[(1, 2), (3, 4)], g[(1, 2, 3)]);
        var (allIm12, allKer12) = AllImKer(a4, g);
        all.UnionWith(allIm12);
        all.UnionWith(allKer12);
        Console.WriteLine($"With A4 subgroups : {all.Count()}");

        // S3 x C2 subgroups and isomorphic are added by hand
        var d12 = Group.Generate(g, g[(1, 2, 3)], g[(1, 2)], g[(4, 5)]);
        var (allIm12b, allKer12b) = AllImKer(d12, g);
        all.UnionWith(allIm12b);
        all.UnionWith(allKer12b);
        Console.WriteLine($"With S3 x C2 subgroups : {all.Count()}");

        // C5 : C4 subgroups and isomorphic are added by hand
        var sdp20 = Group.Generate(g, g[(2, 3, 5, 4)], g[(1, 2, 3, 4, 5)]);
        var (allIm20, allKer20) = AllImKer(sdp20, g);
        all.UnionWith(allIm20);
        all.UnionWith(allKer20);
        Console.WriteLine($"With (C5 : C4) subgroups : {all.Count()}");

        // S4 isomorphic is added by hand
        var s4 = Group.Generate("S4", g, g[(1, 2)], g[(1, 2, 3, 4)]);
        var (allIm24, allKer24) = AllImKer(s4, g);
        all.UnionWith(allIm24);
        all.UnionWith(allKer24);
        Console.WriteLine($"With S4 subgroups : {all.Count()}");

        var allSorted = all.GroupBy(p => p.Count).ToDictionary(a => a.Key, a => a.ToList());
        Console.WriteLine(allSorted.AscendingByKey().Select(p => $"[{p.Key}, {p.Value.Count}]").Glue());
    }
}