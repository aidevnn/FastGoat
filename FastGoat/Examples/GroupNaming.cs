using System.Text.RegularExpressions;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using Group = FastGoat.Structures.Group;

namespace FastGoat.Examples;

public static class GroupNaming
{
    [Flags]
    enum NodeType
    {
        Leaf = 0,
        DirectProduct = 1,
        SemiDirectProduct = 2
    }

    abstract class ITreeElt<T> : IElt<ITreeElt<T>> where T : struct, IElt<T>
    {
        public ConcreteGroup<T>? ContentGroup { get; set; }
        public NodeType ContentType { get; set; }
        public string Name { get; set; } = "C1";
        public int NbDoubleDots => Regex.Count(Name, "x ");
        public int NbStars => Regex.Count(Name, "x:");
        public (GroupType, int, int, string) Infos => (ContentGroup!.GroupType, NbDoubleDots + 100 * NbStars, Name.Length, Name);
        public int Hash => (ContentType, Name).GetHashCode();
        public string NameParenthesis => Name.WithParenthesis();
        public bool Equals(ITreeElt<T>? other) => other?.Name == Name;

        public int CompareTo(ITreeElt<T>? other)
        {
            if (other is null)
                return 1;

            return Infos.CompareTo(other.Infos);
        }

        public override int GetHashCode() => Hash;
    }

    class Leaf<T> : ITreeElt<T> where T : struct, IElt<T>
    {
        public Leaf(AllSubgroups<T> subgroups, DecompType decompType)
        {
            ContentGroup = subgroups.Parent;
            ContentType = NodeType.Leaf;
            if (decompType == DecompType.Abelian)
                Name = Group.AbelianInvariants(ContentGroup).Select(e => e.o).Glue(" x ", "C{0}");
            else if (decompType == DecompType.Extension)
                Name = CommonExtensions(subgroups);
            else if (decompType == DecompType.SimpleNonAbelian)
                Name = SimpleNonAbelians(ContentGroup);
        }

        public override string ToString() => Name;
    }

    class TreeOp<T> : ITreeElt<T> where T : struct, IElt<T>
    {
        public ITreeElt<T> Lhs { get; }
        public ITreeElt<T> Rhs { get; }

        public TreeOp(DecompType t, ITreeElt<T> lhs, ITreeElt<T> rhs, ConcreteGroup<T> g)
        {
            ContentGroup = g;
            ContentType = t == DecompType.DirectProduct ? NodeType.DirectProduct : NodeType.SemiDirectProduct;
            (Lhs, Rhs) = (lhs, rhs);
            if (ContentType == NodeType.DirectProduct)
            {
                if (Lhs.CompareTo(Rhs) == 1)
                    (Lhs, Rhs) = (Rhs, Lhs);

                Name = $"{Lhs.NameParenthesis} x {Rhs.NameParenthesis}";
            }
            else
            {
                Name = $"{Lhs.NameParenthesis} x: {Rhs.NameParenthesis}";
                var name2 = CyclicSdp(ContentGroup, Lhs.ContentGroup!, Rhs.ContentGroup!);
                if (!string.IsNullOrEmpty(name2))
                    Name = name2;
            }
        }

        public override string ToString() => Name;
    }

    [Flags]
    enum DecompType
    {
        Abelian = 0,
        DirectProduct = 1,
        SemiDirectProduct = 2,
        Extension = 3,
        SimpleNonAbelian = 4
    }

    static (AllSubgroups<T> k, AllSubgroups<T> h, DecompType)[] Possibilities<T>(AllSubgroups<T> subgroups) where T : struct, IElt<T>
    {
        var G = subgroups.Parent;
        var tr = subgroups.Restriction(Group.Generate("C1", G, G.Neutral()));
        if (G.GroupType == GroupType.AbelianGroup)
            return new[] { (subgroups, tr, DecompType.Abelian) };

        if (G.GroupType == GroupType.NonAbelianGroup && subgroups.IsSimple())
            return new[] { (subgroups, tr, DecompType.SimpleNonAbelian) };

        var normals = subgroups.Where(sg => sg.IsProperNormal).ToArray();
        var dic = normals.ToDictionary(n => n, n => subgroups.Where(sg => sg.Order == n.Index).ToArray());

        var dirProd = dic.Select(e => (e.Key,
                e.Value.Where(sg => sg.IsNormal && e.Key.Representative.Intersect(sg.Representative).Count() == 1).ToArray()))
            .Where(e => e.Item2.Length != 0)
            .Select(e => (k: e.Key.Representative, h: e.Item2[0].Representative, DecompType.DirectProduct))
            .ToArray();

        var semiDirProd = dic.Select(e => (k: e.Key.Representative, e.Value.Where(sg => !sg.IsNormal).ToArray()))
            .Where(e => e.Item2.Length != 0)
            .Select(e => (e.k, e.Item2.Select(sc => sc.Conjugates.Where(s => s.Intersect(e.k).Count() == 1).ToArray()).ToArray()))
            .Select(e => (e.k, e.Item2.Where(l => l.Length != 0).Select(l => l[0]).ToArray()))
            .SelectMany(e => e.Item2.Select(s => (e.k, h: s, DecompType.SemiDirectProduct)))
            .ToArray();

        var allProds = dirProd.Concat(semiDirProd)
            .Select(e => (subgroups.Restriction(e.k), subgroups.Restriction(e.h), e.Item3))
            .ToArray();

        if (allProds.Length != 0)
            return allProds;

        return new[] { (subgroups, tr, DecompType.Extension) };
    }

    static Homomorphism<T, Automorphism<T>> GroupAction<T>(ConcreteGroup<T> G, ConcreteGroup<T> K, ConcreteGroup<T> H)
        where T : struct, IElt<T>
    {
        if (!G.SuperSetOf(K) || !G.SuperSetOf(H))
            throw new GroupException(GroupExceptionType.NotSubGroup);

        if (K.Intersect(H).Count() > 1)
            throw new();

        var conj = Group.ByConjugate(G);
        if (G.Grid2D(K).Any(e => !K.Contains(conj(e.t1, e.t2))))
            throw new GroupException(GroupExceptionType.NotNormal);

        var map1 = H.ToDictionary(g => g, g => K.ToDictionary(k => k, k => conj(g, k)));
        var autK = Group.AutBase(K);
        var map2 = map1.ToDictionary(e => e.Key, e => new Automorphism<T>(autK, e.Value));
        return new Homomorphism<T, Automorphism<T>>(H, map2);
    }

    static string GroupActionName<T>(Homomorphism<T, Automorphism<T>> act) where T : struct, IElt<T>
    {
        var K = act.Image().First().Domain;
        var H = act.Domain;

        if (act.Image().Distinct().Count() == 1)
            return $"{K.NameParenthesis()} x {H.NameParenthesis()}";

        var gensK = K.GetGenerators().ToArray();
        var gensH = H.GetGenerators().ToArray();

        var (gh, gk) = (gensH[0], gensK[0]);
        var aut = act[gh];
        var dicK = Group.Cycle(K, gk);
        var (m, n, r) = (K.Count(), H.Count(), dicK[aut[gk]]);

        if (n == 2)
        {
            if (m == 3)
                return "S3";
            else if (r == m - 1)
                return $"D{2 * m}";
            else if (int.IsPow2(m))
            {
                if (r == m / 2 - 1)
                    return $"QD{2 * m}";
                else if (r == m / 2 + 1)
                    return $"MM{2 * m}";
            }
        }

        return IntExt.Gcd(m, n * (r - 1)) == 1 ? $"F({m}x:{n}){r}" : $"M({m}x:{n}){r}";
    }

    static string CommonExtensions<T>(AllSubgroups<T> subgroups) where T : struct, IElt<T>
    {
        var G = subgroups.Parent;
        var orders = G.ElementsOrdersList().GroupBy(a => a)
            .ToDictionary(a => a.Key, a => a.Count())
            .AscendingByKey()
            .GlueMap(fmt: "[{0}]:{1}");
        if (G.Count() == 8 && orders == "[1]:1, [2]:1, [4]:6")
            return "Q8";
        if (G.Count() == 16 && orders == "[1]:1, [2]:1, [4]:10, [8]:4")
            return "Q16";
        if (G.Count() == 32 && orders == "[1]:1, [2]:1, [4]:18, [8]:4, [16]:8")
            return "Q32";

        var normals = subgroups.Where(sg => sg.IsProperNormal).ToArray();
        var dic = normals.ToDictionary(n => n, n => subgroups.Where(sg => sg.Order == n.Index).ToArray());
        var candidates = dic.Select(e =>
                (e.Key, e.Value.OrderBy(sc => sc.Representative.GroupType).ThenByDescending(sc => sc.IsMonogenic).ToArray()))
            .OrderByDescending(e => e.Key.GroupType == GroupType.AbelianGroup && e.Item2[0].GroupType == GroupType.AbelianGroup)
            .ThenByDescending(e => e.Key.Order).ThenBy(e => e.Key.GroupType).ToArray();

        if (candidates.Length != 0)
        {
            var c = candidates.First(e => e.Item2[0].Representative.GroupType == GroupType.AbelianGroup);
            var k = BuildName(subgroups.Restriction(c.Key.Representative))[0].Name;
            var h = BuildName(subgroups.Restriction(c.Item2[0].Representative))[0].Name;
            return $"{k.WithParenthesis()} . {h.WithParenthesis()}";
        }

        return $"G[{G.Count()}]";
    }

    static string SimpleNonAbelians<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var og = G.Count();
        if (og == 60)
            return "Alt5";
        if (og == 168)
            return "SL(3,2)";
        if (og == 360)
            return "Alt6";
        if (og == 504)
            return "SL(2,8)";
        if (og == 660)
            return "L2(11)";
        if (og == 1092)
            return "L2(13)";

        return G.Name;
    }

    static string CyclicSdp<T>(ConcreteGroup<T> G, ConcreteGroup<T> K, ConcreteGroup<T> H) where T : struct, IElt<T>
    {
        var gensK = K.GetGenerators().ToArray();
        var gensH = H.GetGenerators().ToArray();

        if (gensK.Length != 1 || gensH.Length != 1)
            return "";

        var act = GroupAction(G, K, H);
        var name = GroupActionName(act);
        return name;
    }

    static ITreeElt<T>[] BuildName<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var subgroups = new AllSubgroups<T>(G);
        return BuildName(subgroups);
    }

    static ITreeElt<T>[] BuildName<T>(AllSubgroups<T> subgroups) where T : struct, IElt<T>
    {
        var all = new List<ITreeElt<T>>();
        var G = subgroups.Parent;
        foreach (var (k, h, t) in Possibilities(subgroups))
        {
            if (h.Count() == 1)
                all.Add(new Leaf<T>(k, t));
            else
                all.AddRange(BuildName(k).Grid2D(BuildName(h)).Select(e => new TreeOp<T>(t, e.t1, e.t2, G)));
        }

        return all.Distinct().Order().ToArray();
    }

    static void ShowNames<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        BuildName(G).Println($"Group:{G.ShortName}");
    }

    public static void Example1()
    {
        ShowNames(FG.Dihedral(4));
        ShowNames(FG.DihedralSdp(5));
        ShowNames(FG.Dihedral(8));
        ShowNames(FG.Symmetric(4));
        ShowNames(Product.Generate(new Cn(5), Group.SemiDirectProd(new Cn(3), new Cn(4))));
        ShowNames(FG.Quaternion(8));
        ShowNames(FG.DiCyclic(6));
        ShowNames(FG.DiCyclic(7));
        ShowNames(FG.DiCyclic(8));
        ShowNames(FG.SemiDihedral(5));
    }

    public static void Example2()
    {
        var allExts24 = FG.AllExtensions((FG.Abelian(12), FG.Abelian(2)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts24)
        {
            var it = BuildName(extInfos.allSubs);
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example3()
    {
        ShowNames(FG.Abelian(5));
        ShowNames(FG.Alternate(5));
        ShowNames(FG.Symmetric(5));
        ShowNames(FG.Abelian(7, 3, 2));
        ShowNames(FG.Abelian(5, 2));
        ShowNames(FG.GLnp(3, 2));
        ShowNames(FG.SL2p(5));

        // Top Down simple non abelian group construction, name dont changes
        ShowNames(FG.L2p(5));
        ShowNames(FG.L2p(7));
        ShowNames(FG.L2p(11));
    }

    public static void Example4()
    {
        var allExts20 = FG.AllExtensions(FG.AllAbelianGroupsOfOrder(4).Select(e => (FG.Abelian(5), e)).ToArray())
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        var allExts40 = FG.AllExtensions(FG.AllAbelianGroupsOfOrder(4).Select(e => (FG.Abelian(2, 5), e)).ToArray())
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts20)
        {
            var it = BuildName(extInfos.allSubs);
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }

        Console.WriteLine();

        foreach (var extInfos in allExts40)
        {
            var it = BuildName(extInfos.allSubs);
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example5()
    {
        var allExts = FG.AllExtensions((FG.Abelian(4, 4), FG.Abelian(2)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = BuildName(extInfos.allSubs);
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example6()
    {
        var allExts = FG.AllExtensions(
                (FG.Abelian(2, 8), FG.Abelian(2)),
                (FG.Abelian(4, 4), FG.Abelian(2)),
                (FG.Abelian(2, 4), FG.Abelian(4)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = BuildName(extInfos.allSubs);
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }

    public static void Example7()
    {
        var allExts = FG.AllExtensions((FG.SL2p(3), FG.Abelian(2)))
            .OrderBy(e => e.ext.GroupType)
            .ThenByDescending(e => e.ext.ElementsOrders.Values.Max())
            .ThenBy(e => ((int, int, int))e.allSubs.Infos).ToList();

        foreach (var extInfos in allExts)
        {
            var it = BuildName(extInfos.allSubs);
            extInfos.ext.Name = it.First().Name;
            CocyclesDFS.DisplayInfosGroups([(extInfos.ext, ((int, int, int))extInfos.allSubs.Infos)], naming: false);
            it.Println("Group Names");
        }
    }
}