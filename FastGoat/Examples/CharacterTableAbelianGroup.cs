using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;

namespace FastGoat.Examples;

public static class CharacterTableAbelianGroup
{
    public struct StringElt : IElt<StringElt>
    {
        public string Content { get; }

        public StringElt()
        {
            Content = "";
        }

        public StringElt(object o)
        {
            Content = o.ToString() ?? string.Empty;
        }

        public bool Equals(StringElt other) => Content.Equals(other.Content);

        public int CompareTo(StringElt other) => String.CompareOrdinal(Content, other.Content);

        public int Hash => Content.GetHashCode();
        public override string ToString() => Content;
        public override int GetHashCode() => Hash;
    }

    static void CheckProperties<T>(ConcreteGroup<T> g, Dictionary<T, KMatrix<EPoly<Rational>>> table) where T : struct, IElt<T>
    {
        var e0 = table.Values.First().ToArray().First();
        var n = g.Count();
        var rg = n.Range();
        var keys = table.Keys.ToArray();
        var allCombs = rg.SelectMany(i => rg.Where(j => j > i).Select(j => (keys[i], g.Invert(keys[j])))).ToArray();
        Console.WriteLine("Sum[Ꭓi(g)Ꭓj(g^−1 )]= 0  : {0}", allCombs.All(e => (table[e.Item1].T * table[e.Item2]).IsZero()));
        Console.WriteLine("Sum[Ꭓi(g)Ꭓi(g^−1 )]=|G| : {0}", rg.All(i => (table[keys[i]].T * table[g.Invert(keys[i])])[0, 0].Equals(e0.One * n)));
    }

    static void TableCn(int n)
    {
        var nth = new NthRootQ(n);
        var w = nth.PrimitivesRoots().First();
        EPoly<Rational> Chi(int r, ZnInt a) => w.Pow(r * a.K);

        var table = new StringElt[n + 1, n + 1];
        var mat = new KMatrix<EPoly<Rational>>(w.Zero, n, n);
        var cn = new Cn(n);
        for (int i = 0; i < n; i++)
        {
            table[0, i + 1] = new($"Cl({cn[i]})");
            table[i + 1, 0] = new($"Ꭓ{i + 1}");
            for (int j = 0; j < n; j++)
            {
                var c = mat.Coefs[i, j] = Chi(i, cn[j]);
                table[i + 1, j + 1] = new(c);
            }
        }

        Ring.MatrixDisplayForm = Ring.MatrixDisplay.Table;
        DisplayGroup.Head(cn);
        Console.WriteLine($"With {w.F} = 0");
        Ring.DisplayMatrix(table);

        var tableDic = n.Range().ToDictionary(i => cn[i], i => mat.GetCol(i));
        CheckProperties(cn, tableDic);
        Console.WriteLine();
    }

    static void TableCmCn(int m, int n)
    {
        var nth = new NthRootQ(m * n);
        var w = nth.PrimitivesRoots().First();
        var wm = w.Pow(n);
        var wn = w.Pow(m);

        EPoly<Rational> Chi(int r, Ep<ZnInt> a)
        {
            var (rn, rm) = Int32.DivRem(r, m);
            return (wm.Pow(a[0].K * rm) * wn.Pow(a[1].K * rn));
        }

        var table = new StringElt[m * n + 1, m * n + 1];
        var cmcn = FG.Abelian(m, n);
        var mat = new KMatrix<EPoly<Rational>>(w.Zero, m * n, m * n);
        for (int ni = 0; ni < n; ni++)
        {
            for (int mi = 0; mi < m; mi++)
            {
                var ri = mi + m * ni;

                table[0, ri + 1] = new($"Cl{cmcn[mi, ni]}");
                table[ri + 1, 0] = new($"Ꭓ{ri + 1}");

                for (int nj = 0; nj < n; nj++)
                {
                    for (int mj = 0; mj < m; mj++)
                    {
                        var rj = mj + m * nj;
                        var c = mat.Coefs[ri, rj] = Chi(ri, cmcn[mj, nj]);
                        table[ri + 1, rj + 1] = new(c);
                    }
                }
            }
        }

        Ring.MatrixDisplayForm = Ring.MatrixDisplay.Table;
        DisplayGroup.Head(cmcn);
        Console.WriteLine($"With {w.F} = 0");
        Ring.DisplayMatrix(table);

        var tableDic = m.Range().Grid2D(n.Range()).ToDictionary(i => cmcn[i.t1, i.t2], i => mat.GetCol(i.t1 + m * i.t2));
        CheckProperties(cmcn, tableDic);
        Console.WriteLine();
    }

    public static void ConjugacyClasses()
    {
        Group.DisplayConjugacyClasses(FG.Abelian(2, 2));
        Group.DisplayConjugacyClasses(FG.Abelian(2, 3));
        Group.DisplayConjugacyClasses(FG.Abelian(3, 3));
        Group.DisplayConjugacyClasses(FG.Abelian(8));
        Group.DisplayConjugacyClasses(FG.Abelian(3, 4));
        Group.DisplayConjugacyClasses(FG.Abelian(15));
        Group.DisplayConjugacyClasses(FG.Abelian(3, 5));
    }

    public static void ExampleTableCn()
    {
        TableCn(2);
        TableCn(3);
        TableCn(4);
        TableCn(5);
        TableCn(6);
        TableCn(8);
        TableCn(12);

        TableCmCn(2, 2);
        TableCmCn(2, 3);
        TableCmCn(2, 4);
        TableCmCn(3, 4);
    }

    public static void ExamplesCharactersTableAbelianGroups()
    {
        FG.CharactersTable(FG.Dihedral(8)).DisplayCells();
        FG.CharactersTable(FG.Abelian(3)).DisplayCells();
        FG.CharactersTable(FG.Abelian(2, 2)).DisplayCells();
        FG.CharactersTable(FG.Abelian(2, 3)).DisplayCells();
        FG.CharactersTable(FG.Abelian(6)).DisplayCells();
        FG.CharactersTable(FG.Abelian(4)).DisplayCells();
        FG.CharactersTable(FG.Abelian(2, 4)).DisplayCells();
        FG.CharactersTable(FG.Abelian(2, 2, 2)).DisplayCells();
        FG.CharactersTable(FG.Abelian(2, 2, 3)).DisplayCells();
    }

    public static void ExamplesLiftDerivedGroup()
    {
        FG.CharactersTable(FG.Symmetric(3)).DisplayCells();
        FG.CharactersTable(FG.Dihedral(4)).DisplayCells();
        FG.CharactersTable(FG.Alternate(4)).DisplayCells();
        FG.CharactersTable(FG.Quaternion(8)).DisplayCells();
        FG.CharactersTable(FG.Dihedral(5)).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(new Cn(5), new Cn(4))).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(new Cn(7), new Cn(6))).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(new Cn(11), new Cn(10))).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(new Cn(13), new Cn(12))).DisplayCells();
    }

    public static void ExamplesTwoMissingCharacters()
    {
        FG.CharactersTable(Group.SemiDirectProd(new Cn(3), new Cn(4))).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(new Cn(4), new Cn(4))).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(new Cn(7), new Cn(3))).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(new Cn(9), new Cn(3))).DisplayCells();
        FG.CharactersTable(FG.Dihedral(6)).DisplayCells();
        FG.CharactersTable(FG.DiCyclic(3)).DisplayCells();
        FG.CharactersTable(FG.SemiDihedral(4)).DisplayCells();
        FG.CharactersTable(FG.Symmetric(4)).DisplayCells();
        
        FG.CharactersTable(Group.SemiDirectProd(FG.Abelian(3, 3), new Cn(4))).DisplayCells();
        // FG.CharactersTable(new WordGroup("a3, b3, c4, ab = ba, cac−1 = b, cbc−1 = a2")).DisplayCells();
        // FG.CharactersTable(Group.SemiDirectProd(new Cn(3), new Cn(8))).DisplayCells(); // wont work
    }

    public static void ExamplesPQgroups()
    {
        for (int i = 2; i <= 16; i++)
            FG.CharactersTable(FG.Dihedral(i)).DisplayCells();
        
        var pqGroups = 31.Range(2).SelectMany(i => FG.Frobenius(i)).ToArray();
        foreach (var g in pqGroups)
        {
            FG.CharactersTable(g).DisplayCells();
        }
    }
}