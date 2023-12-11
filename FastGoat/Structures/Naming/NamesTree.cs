using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Naming;

public static class NamesTree
{
    static Leaf[] CommonNames(ConcreteGroup<TableElt> G)
    {
        var og = G.Count();
        var orders = G.ElementsOrdersList().GroupBy(a => a)
            .ToDictionary(a => a.Key, a => a.Count())
            .AscendingByKey()
            .GlueMap(fmt: "[{0}]:{1}");

        if (og == 12 && orders == "[1]:1, [2]:3, [3]:8")
            return [new Leaf(G, "A4")];
        if (og == 24 && orders == "[1]:1, [2]:9, [3]:8, [4]:6")
            return [new Leaf(G, "S4")];
        if (og == 120 && orders == "[1]:1, [2]:25, [3]:20, [4]:30, [5]:24, [6]:20")
            return [new Leaf(G, "S5")];
        if (og == 720 && orders == "[1]:1, [2]:75, [3]:80, [4]:180, [5]:144, [6]:240")
            return [new Leaf(G, "S6")];
        
        if (og == 8 && orders == "[1]:1, [2]:1, [4]:6")
            return [new Leaf(G, "Q8")];
        if (og == 16 && orders == "[1]:1, [2]:1, [4]:10, [8]:4")
            return [new Leaf(G, "Q16")];
        if (og == 32 && orders == "[1]:1, [2]:1, [4]:18, [8]:4, [16]:8")
            return [new Leaf(G, "Q32")];
        
        if (og == 24 && orders == "[1]:1, [2]:1, [3]:8, [4]:6, [6]:8")
            return [new Leaf(G, "SL(2,3)")];
        if (og == 48 && orders == "[1]:1, [2]:13, [3]:8, [4]:6, [6]:8, [8]:12")
            return [new Leaf(G, "GL(2,3)")];
        if (og == 120 && orders == "[1]:1, [2]:1, [3]:20, [4]:30, [5]:24, [6]:20, [10]:24")
            return [new Leaf(G, "SL(2,5)")];
        if (og == 480 && orders == "[1]:1, [2]:31, [3]:20, [4]:152, [5]:24, [6]:20, [8]:40, [10]:24, [12]:40, [20]:48, [24]:80")
            return [new Leaf(G, "GL(2,5)")];
        
        if (og == 12 && orders == "[1]:1, [2]:1, [3]:2, [4]:6, [6]:2")
            return [new Leaf(G, "Dic3")];
        if (og == 16 && orders == "[1]:1, [2]:1, [4]:10, [8]:4")
            return [new Leaf(G, "Dic4")];
        if (og == 20 && orders == "[1]:1, [2]:1, [4]:10, [5]:4, [10]:4")
            return [new Leaf(G, "Dic5")];
        if (og == 24 && orders == "[1]:1, [2]:1, [3]:2, [4]:14, [6]:2, [12]:4")
            return [new Leaf(G, "Dic6")];
        if (og == 28 && orders == "[1]:1, [2]:1, [4]:14, [7]:6, [14]:6")
            return [new Leaf(G, "Dic7")];
        if (og == 36 && orders == "[1]:1, [2]:1, [3]:2, [4]:18, [6]:2, [9]:6, [18]:6")
            return [new Leaf(G, "Dic9")];
        if (og == 40 && orders == "[1]:1, [2]:1, [4]:22, [5]:4, [10]:4, [20]:8")
            return [new Leaf(G, "Dic10")];
        if (og == 44 && orders == "[1]:1, [2]:1, [4]:22, [11]:10, [22]:10")
            return [new Leaf(G, "Dic11")];
        if (og == 48 && orders == "[1]:1, [2]:1, [3]:2, [4]:26, [6]:2, [8]:4, [12]:4, [24]:8")
            return [new Leaf(G, "Dic12")];
        if (og == 52 && orders == "[1]:1, [2]:1, [4]:26, [13]:12, [26]:12")
            return [new Leaf(G, "Dic13")];
        if (og == 56 && orders == "[1]:1, [2]:1, [4]:30, [7]:6, [14]:6, [28]:12")
            return [new Leaf(G, "Dic14")];
        if (og == 60 && orders == "[1]:1, [2]:1, [3]:2, [4]:30, [5]:4, [6]:2, [10]:4, [15]:8, [30]:8")
            return [new Leaf(G, "Dic15")];
        if (og == 68 && orders == "[1]:1, [2]:1, [4]:34, [17]:16, [34]:16")
            return [new Leaf(G, "Dic17")];
        if (og == 72 && orders == "[1]:1, [2]:1, [3]:2, [4]:38, [6]:2, [9]:6, [12]:4, [18]:6, [36]:12")
            return [new Leaf(G, "Dic18")];
        if (og == 76 && orders == "[1]:1, [2]:1, [4]:38, [19]:18, [38]:18")
            return [new Leaf(G, "Dic19")];
        if (og == 80 && orders == "[1]:1, [2]:1, [4]:42, [5]:4, [8]:4, [10]:4, [20]:8, [40]:16")
            return [new Leaf(G, "Dic20")];
        if (og == 84 && orders == "[1]:1, [2]:1, [3]:2, [4]:42, [6]:2, [7]:6, [14]:6, [21]:12, [42]:12")
            return [new Leaf(G, "Dic21")];
        if (og == 88 && orders == "[1]:1, [2]:1, [4]:46, [11]:10, [22]:10, [44]:20")
            return [new Leaf(G, "Dic22")];
        if (og == 92 && orders == "[1]:1, [2]:1, [4]:46, [23]:22, [46]:22")
            return [new Leaf(G, "Dic23")];
        if (og == 96 && orders == "[1]:1, [2]:1, [3]:2, [4]:50, [6]:2, [8]:4, [12]:4, [16]:8, [24]:8, [48]:16")
            return [new Leaf(G, "Dic24")];

        return Array.Empty<Leaf>();
    }

    static (AllSubgroups<TableElt> k, AllSubgroups<TableElt> h, ANameElt.DecompType)[] AllOps(AllSubgroups<TableElt> subgroups)
    {
        var G = subgroups.Parent;
        var tr = subgroups.Restriction(Group.Generate("C1", G, G.Neutral()));
        if (G.GroupType == GroupType.AbelianGroup)
            return new[] { (subgroups, tr, ANameElt.DecompType.Abelian) };

        if (G.GroupType == GroupType.NonAbelianGroup && subgroups.IsSimple())
            return new[] { (subgroups, tr, ANameElt.DecompType.SimpleNonAbelian) };

        var normals = subgroups.Where(sg => sg.IsProperNormal).ToArray();
        var dic = normals.ToDictionary(n => n, n => subgroups.Where(sg => sg.Order == n.Index).ToArray());

        var dirProd = dic.Select(e => (e.Key,
                e.Value.Where(sg => sg.IsNormal && e.Key.Representative.Intersect(sg.Representative).Count() == 1).ToArray()))
            .Where(e => e.Item2.Length != 0)
            .Select(e => (k: e.Key.Representative, h: e.Item2[0].Representative, ANameElt.DecompType.DirectProduct))
            .ToArray();

        var semiDirProd = dic.Select(e => (k: e.Key.Representative, e.Value.Where(sg => !sg.IsNormal).ToArray()))
            .Where(e => e.Item2.Length != 0)
            .Select(e => (e.k, e.Item2.Select(sc => sc.Conjugates.Where(s => s.Intersect(e.k).Count() == 1).ToArray()).ToArray()))
            .Select(e => (e.k, e.Item2.Where(l => l.Length != 0).Select(l => l[0]).ToArray()))
            .SelectMany(e => e.Item2.Select(s => (e.k, h: s, ANameElt.DecompType.SemiDirectProduct)))
            .ToArray();

        var allProds = dirProd.Concat(semiDirProd)
            .Select(e => (subgroups.Restriction(e.k), subgroups.Restriction(e.h), e.Item3))
            .ToArray();

        var usedNormals = dirProd.SelectMany(e => new[] { e.k, e.h }).Concat(semiDirProd.Select(e => e.k))
            .Select(sg => subgroups.First(sc => sc.Representative.SetEquals(sg))).ToArray();
        var remNormals = normals.Except(usedNormals).ToArray();
        var extOps = remNormals.Select(e => (e, G.Over(e.Representative).ToTable()))
            .Select(e => (subgroups.Restriction(e.e.Representative), new AllSubgroups<TableElt>(e.Item2),
                ANameElt.DecompType.Extension))
            .ToArray();

        return [..allProds, ..extOps];
    }

    public static ANameElt[] BuildName(AllSubgroups<TableElt> subgroups)
    {
        var all = new List<ANameElt>();
        var G = subgroups.Parent;
        if (G.Count() > 3000)
            throw new GroupException(GroupExceptionType.GroupDef);

        foreach (var (k, h, t) in AllOps(subgroups))
        {
            if (t == ANameElt.DecompType.Abelian || t == ANameElt.DecompType.SimpleNonAbelian)
                all.Add(new Leaf(k, t));
            else if (t == ANameElt.DecompType.SemiDirectProduct)
                all.AddRange(BuildName(k).Grid2D(BuildName(h)).Select(e => new SemiDirectProductOp(e.t1, e.t2, G)).Distinct()
                    .Order().Take(1));
            else if (t == ANameElt.DecompType.DirectProduct)
                all.AddRange(BuildName(k).Grid2D(BuildName(h)).Select(e => new DirectProductOp(e.t1, e.t2, G)).Distinct().Order()
                    .Take(1));
            else
                all.AddRange(BuildName(k).Grid2D(BuildName(h)).Select(e => new ExtensionOp(e.t1, e.t2, G)).Distinct().Order()
                    .Take(1));
        }

        return all.Concat(CommonNames(G)).Distinct().Order().ToArray();
    }

    public static ANameElt[] BuildName<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var subGroups = new AllSubgroups<TableElt>(G.ToTable());
        return BuildName(subGroups);
    }
}