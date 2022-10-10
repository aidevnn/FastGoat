using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class AllSubGroups
{
    static (HashSet<HashSet<T>> im, HashSet<HashSet<T>> ker) AllImKer<T>(ConcreteGroup<T> g0, ConcreteGroup<T> g1)
        where T : struct, IElt<T>
    {
        HashSet<HashSet<T>> all = new HashSet<HashSet<T>>(new SeqEquality<T>());
        HashSet<HashSet<T>> allIm = new HashSet<HashSet<T>>(new SeqEquality<T>());
        HashSet<HashSet<T>> allKer = new HashSet<HashSet<T>>(new SeqEquality<T>());
        var hom = Group.AllHomomorphisms(g0, g1);
        foreach (var m in hom)
        {
            var im = m.Values.ToHashSet();
            var ker = m.Where(p => p.Value.Equals(g0.Neutral())).Select(p => p.Key).ToHashSet();
            allIm.Add(im);
            allKer.Add(ker);
        }

        return (allIm, allKer);
    }

    static void AllNormalSubGroups<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var allHoms = Group.AllHomomorphisms(g, g);
        var allKer = allHoms.Select(p => p.Where(kp => kp.Value.Equals(g.Neutral())).Select(kp => kp.Key).ToHashSet())
            .ToHashSet(new SeqEquality<T>());
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
        var g = Group.Generate("A4", s4[(1, 3), (2, 4)], s4[(1, 2, 3)]);
        AllNormalSubGroups(g);
    }

    public static void A5Simple()
    {
        var s5 = new Sn(5);
        var g = Group.Generate("A5", s5[(1, 3, 5)], s5[(1, 2, 3, 4, 5)]);
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
        var all = allIm.Union(allKer).ToHashSet(new SeqEquality<Ep2<ZnInt, ZnInt>>());
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void AbelianSubGroupsAnother()
    {
        var g = Product.Generate(new Cn(6), new Cn(9));
        var (allIm, allKer) = AllImKer(g, g);
        var all = allIm.Union(allKer).ToHashSet(new SeqEquality<Ep2<ZnInt, ZnInt>>());
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void AbelianSubGroupsEncore()
    {
        var g = Product.Generate(new Cn(2), new Cn(4), new Cn(8));
        var (allIm, allKer) = AllImKer(g, g);
        var all = allIm.Union(allKer).ToHashSet(new SeqEquality<Ep3<ZnInt, ZnInt, ZnInt>>());
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void Symm3SubGroups()
    {
        var g = Group.Create(new Sn(3));
        var all = new HashSet<HashSet<Perm>>(new SeqEquality<Perm>());
        all.UnionWith(g.Select(e => Group.Generate(g, e).ToHashSet()));
        Console.WriteLine($"Monogenics SubGroups from elements : {all.Count()}");

        var (_, allKer) = AllImKer(g, g);
        all.UnionWith(allKer);
        Console.WriteLine($"{g} AllSubGroups : {all.Count()}");
    }

    public static void Symm4SubGroups()
    {
        var g = Group.Create(new Sn(4));
        var all = new HashSet<HashSet<Perm>>(new SeqEquality<Perm>());
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
}