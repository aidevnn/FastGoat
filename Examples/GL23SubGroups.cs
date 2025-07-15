using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace Examples;

public static class GL23SubGroups
{
    // Problem 6.23 (Project—Subgroups of GL2(3)) page 137, H.E. Rose
    static GL23SubGroups()
    {
        var gl23 = new GL(2, 3);
        var r0 = gl23[2, 1, 0, 1];
        var r1 = gl23[1, 0, 1, 1];
        GL23mat = Group.Generate("GL2(3)", gl23, r0, r1);
        GL23wg = new WordGroup("a8, b2, c3, bab = a3, bcb = c2, c2a2c = ab, c2abc = aba2");
    }

    private static WordGroup GL23wg { get; }
    private static ConcreteGroup<Mat> GL23mat { get; }

    public static void ShowAll()
    {
        var all = new HashSet<ConcreteGroup<Mat>>(new GroupSetEquality<Mat>());
        all.UnionWith(CyclicSubGroups());
        all.UnionWith(KleinSubGroups());
        all.UnionWith(Symm3SubGroups());
        all.Add(Quaternion());
        all.UnionWith(Dihedral8SubGroups());
        all.UnionWith(Order12Subgroups());
        all.UnionWith(SylowPSubgroups());
        all.Add(SL23());
        all.Add(GL23mat);

        foreach (var kp0 in all.GroupBy(g => g.Count()).OrderBy(e => e.Key))
        {
            Console.WriteLine($"Order : {kp0.Key} Count : {kp0.Count()}");
            foreach (var kp1 in kp0.GroupBy(g => g.GroupType))
            {
                Console.WriteLine($"  {kp1.Key}");
                Console.WriteLine($"    {kp1.OrderBy(g => g.Name).Glue(", ")}");
            }
        }

        Console.WriteLine();
        Console.WriteLine($"Total : {all.Count}");
    }

    private static ConcreteGroup<Mat> Quaternion()
    {
        var z1 = Group.Zentrum(GL23mat);
        var ord4H1 = GL23mat.Where(e => GL23mat.ElementsOrders[e] == 4).ToArray();
        var h1 = Group.Generate("Q8", GL23mat, z1.Concat(ord4H1).ToArray());
        return h1;
    }

    private static ConcreteGroup<Mat> SL23()
    {
        var gl23 = new GL(2, 3);
        var sl23 = GL23mat.Where(e => gl23.Det(e) == 1).ToArray();
        var h1 = Group.Generate("SL2(3)", GL23mat, sl23.ToArray());
        return h1;
    }

    private static List<ConcreteGroup<Mat>> CyclicSubGroups()
    {
        var allCyclics = GL23mat.Select(e => Group.Generate(GL23mat, e)).ToHashSet(new GroupSetEquality<Mat>());
        var groups = allCyclics.GroupBy(a => a.Count()).ToDictionary(e => e.Key, e => e.ToList());
        foreach (var (o, g) in groups)
        {
            var i = 1;
            g.ForEach(g0 => g0.Name = $"C{o}[{i++,2:00}]");
        }

        return groups.SelectMany(g => g.Value).ToList();
    }

    private static List<ConcreteGroup<Mat>> KleinSubGroups()
    {
        var g0 = Product.Generate(new Cn(2), new Cn(2));
        return Group.IsomorphicsSubgroupsAll(GL23mat, g0, "Klein");
    }

    private static List<ConcreteGroup<Mat>> Symm3SubGroups()
    {
        var g0 = new Symm(3);
        return Group.IsomorphicsSubgroupsAll(GL23mat, g0, "Symm3");
    }

    private static List<ConcreteGroup<Mat>> SylowPSubgroups()
    {
        var h1 = Quaternion();
        return CyclicSubGroups().Where(e => e.Count() == 8)
            .Select((s, i) => Group.Generate($"2Sylows[{i}]", GL23mat, h1.Concat(s).ToArray())).ToList();
    }

    private static List<ConcreteGroup<Mat>> Dihedral8SubGroups()
    {
        var d8 = Group.SemiDirectProd(new Cn(4), new Cn(2));
        return Group.IsomorphicsSubgroupsAll(GL23mat, d8);
    }

    private static List<ConcreteGroup<Mat>> Order12Subgroups()
    {
        var h12u = Group.Generate("H12", GL23mat, GL23mat[2, 1, 0, 1], GL23mat[1, 2, 0, 1], GL23mat[1, 1, 0, 2]);
        return Group.IsomorphicsSubgroupsAll(GL23mat, h12u);
    }

    public static void ShowGL32Elements()
    {
        DisplayGroup.HeadElements(GL23mat);
        DisplayGroup.HeadElements(GL23wg);
        DisplayGroup.AreIsomorphics(GL23mat, GL23wg);
    }

    public static void ShowCyclicSubGroups()
    {
        foreach (var subgroup in CyclicSubGroups())
        {
            DisplayGroup.HeadOrders(subgroup);
        }
    }

    public static void ShowSymm3SubGroups()
    {
        var all = Symm3SubGroups();
        foreach (var subgroup in all)
        {
            DisplayGroup.HeadElements(subgroup);
        }

        var g1 = all[0];
        var conj1 = Group.SubGroupsConjugates(GL23mat, g1);

        Console.WriteLine($"Nb Inner Conjugates for {g1} : {conj1.Count}");
        Console.WriteLine($"Nb Isomorphics for {g1} : {all.Count}");
    }

    public static void ShowQuaternion()
    {
        var h1 = Quaternion();
        DisplayGroup.HeadElements(h1);

        var q = new WordGroup("Q8", "a4, a2 = b2, ba = a3b");
        DisplayGroup.HeadElements(q);
        DisplayGroup.AreIsomorphics(h1, q);
        Console.WriteLine();

        DisplayGroup.Head(GL23mat.Over(h1));
    }

    public static void ShowSylowSubGroups()
    {
        foreach (var subgroup in SylowPSubgroups())
        {
            DisplayGroup.HeadOrders(subgroup);
        }
    }

    public static void ShowOrder12SubGroups()
    {
        foreach (var subgroup in Order12Subgroups())
        {
            DisplayGroup.HeadElements(subgroup);
        }
    }

    public static void ShowKleinSubGroups()
    {
        foreach (var subgroup in KleinSubGroups())
        {
            DisplayGroup.HeadElements(subgroup);
        }
    }

    public static void ShowD8SubGroups()
    {
        foreach (var subgroup in Dihedral8SubGroups())
        {
            DisplayGroup.HeadElements(subgroup);
        }
    }
}