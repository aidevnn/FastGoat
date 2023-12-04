using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class GroupNaming
{
    [Flags]
    enum TreeContent
    {
        Leaf = 0,
        DirectProduct = 1,
        SemiDirectProduct = 2
    }

    abstract class ITreeElt : IElt<ITreeElt>
    {
        public TreeContent Content { get; set; }
        public string Name { get; set; } = "()";
        public int Hash => (Content, Name).GetHashCode();
        public string NameParenthesis => Name.WithParenthesis();
        public bool Equals(ITreeElt? other) => other?.Name == Name;
        public int CompareTo(ITreeElt? other)
        {
            if (other is null)
                return 1;

            return String.Compare(Name, other.Name, StringComparison.CurrentCulture);
        }

        public override int GetHashCode() => Hash;
    }

    class Leaf<T> : ITreeElt where T : struct, IElt<T>
    {
        public Leaf(ConcreteGroup<T> g)
        {
            G = g;
            Content = TreeContent.Leaf;
            if (g.GroupType == GroupType.AbelianGroup)
                Name = Group.AbelianInvariants(G).Select(e => e.o).Glue(" x ", "C{0}");
            else
                Name = $"G[{g.Count()}]";
        }

        public ConcreteGroup<T> G { get; }
        public override string ToString() => Name;
    }

    class TreeOp<T> : ITreeElt where T : struct, IElt<T>
    {
        public ITreeElt Lhs { get; }
        public ITreeElt Rhs { get; }

        public TreeOp(bool isDirectProd, ITreeElt lhs, ITreeElt rhs)
        {
            var treeContent = isDirectProd ? TreeContent.DirectProduct : TreeContent.SemiDirectProduct;
            (Content, Lhs, Rhs) = (treeContent, lhs, rhs);
            if (Content == TreeContent.DirectProduct)
            {
                if (Lhs.CompareTo(Rhs) == -1)
                    (Lhs, Rhs) = (Rhs, Lhs);
                
                Name = $"{Lhs.NameParenthesis} x {Rhs.NameParenthesis}";
            }
            else 
                Name = $"{Lhs.NameParenthesis} x: {Rhs.NameParenthesis}";
        }

        public override string ToString() => Name;
    }

    [Flags]
    enum DecompType
    {
        DirectProduct = 0,
        SemiDirectProduct = 1,
        Extension = 2
    }

    static (ConcreteGroup<T> k, ConcreteGroup<T> h, DecompType)[] Possibilities<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var tr = Group.Generate("()", G, G.Neutral());
        if (G.GroupType == GroupType.AbelianGroup)
            return new[] { (G, tr, DecompType.DirectProduct) };

        var subs = Group.AllSubGroups(G);
        var subgroups = new AllSubgroups<T>(subs);
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

        var allProds = dirProd.Concat(semiDirProd).ToArray();
        if (allProds.Length != 0)
            return allProds;

        return new[] { (G, tr, DecompType.Extension) };
    }

    static ITreeElt[] BuildName<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        var all = new List<ITreeElt>();
        foreach (var (k, h, t) in Possibilities(G))
        {
            if (h.Count() == 1)
                all.Add(new Leaf<T>(k));
            else
                all.AddRange(BuildName(k).Grid2D(BuildName(h)).Select(e => new TreeOp<T>(t == DecompType.DirectProduct, e.t1, e.t2)));
        }

        return all.Distinct().ToArray();
    }

    static void ShowNames<T>(ConcreteGroup<T> G) where T : struct, IElt<T>
    {
        BuildName(G).Println($"Group:{G.Name}");
    }

    public static void Example1()
    {
        ShowNames(FG.Dihedral(8));
        ShowNames(FG.DihedralSdp(5));
        ShowNames(FG.Symmetric(4));
        ShowNames(Product.Generate(new Cn(5), Group.SemiDirectProd(new Cn(3), new Cn(4))));
        ShowNames(FG.Quaternion(8));
        ShowNames(FG.DiCyclic(6));
        ShowNames(FG.DiCyclic(7));
        ShowNames(FG.DiCyclic(8));
        ShowNames(FG.SemiDihedral(5));
    }
}