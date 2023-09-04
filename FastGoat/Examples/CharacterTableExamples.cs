using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class CharacterTableExamples
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
        Console.WriteLine("Sum[Ꭓi(g)Ꭓi(g^−1 )]=|G| : {0}",
            rg.All(i => (table[keys[i]].T * table[g.Invert(keys[i])])[0, 0].Equals(e0.One * n)));
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
        FG.CharactersTable(FG.Dihedral(8)).DisplayCells();
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
        FG.CharactersTable(Group.SemiDirectProd(new Cn(3), new Cn(8))).DisplayCells();
        FG.CharactersTable(Group.SemiDirectProd(FG.Abelian(3, 3), new Cn(4))).DisplayCells();

        FG.CharactersTable(FG.Dihedral(6)).DisplayCells();
        FG.CharactersTable(FG.DiCyclic(3)).DisplayCells();
        FG.CharactersTable(FG.SemiDihedral(4)).DisplayCells();
    }

    public static void ExamplesPQgroups()
    {
        for (int i = 2; i <= 16; i++)
            FG.CharactersTable2(FG.DihedralSdp(i)).DisplayCells();

        var pqGroups = 31.Range(2).SelectMany(i => FG.FrobeniusSdp(i)).ToArray();
        foreach (var g in pqGroups)
        {
            FG.CharactersTable2(g).DisplayCells();
        }
    }

    public static void ExamplesPermGroups()
    {
        GlobalStopWatch.Restart();
        for (int n = 3; n < 7; n++)
        {
            GlobalStopWatch.AddLap();
            var ctAn = FG.CharactersTable2Slow(FG.Alternate(n));
        
            if (n == 6)
                ctAn.InductionFromSubGroups(2);
            ctAn.DisplayCells();
            GlobalStopWatch.Show($"{ctAn.Gr}");
        }

        for (int n = 3; n < 7; n++)
        {
            GlobalStopWatch.AddLap();
            var ctSn = FG.CharactersTable2Slow(FG.Symmetric(n));
            ctSn.DisplayCells();
            GlobalStopWatch.Show($"{ctSn.Gr}");
        }
    }
    /*
        |Alt5| = 60
        Type        NonAbelianGroup
        BaseGroup   S5
        SuperGroup  |Symm5| = 120

        [Class     1  2  3            5a            5b]
        [ Size     1 15 20            12            12]
        [                                             ]
        [  X.1     1  1  1             1             1]
        [  X.2     3 -1  0   -ξ5³ + -ξ5² ξ5³ + ξ5² + 1]
        [  X.3     3 -1  0 ξ5³ + ξ5² + 1   -ξ5³ + -ξ5²]
        [  X.4     4  0  1            -1            -1]
        [  X.5     5  1 -1             0             0]
        All i,                 Sum[g](Xi(g)Xi(g^−1))= |G|      : True
        All i <> j,            Sum[g](Xi(g)Xj(g^−1))=  0       : True
        All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|  : True
        All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0       : True

        |Symm5| = 120
        Type        NonAbelianGroup
        BaseGroup   S5

        [Class     1 2a 2b  3  4  5  6]
        [ Size     1 10 15 20 30 24 20]
        [                             ]
        [  X.1     1  1  1  1  1  1  1]
        [  X.2     1 -1  1  1 -1  1 -1]
        [  X.3     4  2  0  1  0 -1 -1]
        [  X.4     4 -2  0  1  0 -1  1]
        [  X.5     5  1  1 -1 -1  0  1]
        [  X.6     5 -1  1 -1  1  0 -1]
        [  X.7     6  0 -2  0  0  1  0]
        All i,                 Sum[g](Xi(g)Xi(g^−1))= |G|      : True
        All i <> j,            Sum[g](Xi(g)Xj(g^−1))=  0       : True
        All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|  : True
        All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0       : True

        |Alt6| = 360
        Type        NonAbelianGroup
        BaseGroup   S6
        SuperGroup  |Symm6| = 720

        [Class      1  2 3a 3b  4            5a            5b]
        [ Size      1 45 40 40 90            72            72]
        [                                                    ]
        [  X.1      1  1  1  1  1             1             1]
        [  X.2      5  1  2 -1 -1             0             0]
        [  X.3      5  1 -1  2 -1             0             0]
        [  X.4      8  0 -1 -1  0   -ξ5³ + -ξ5² ξ5³ + ξ5² + 1]
        [  X.5      8  0 -1 -1  0 ξ5³ + ξ5² + 1   -ξ5³ + -ξ5²]
        [  X.6      9  1  0  0  1            -1            -1]
        [  X.7     10 -2  1  1  0             0             0]
        All i,                 Sum[g](Xi(g)Xi(g^−1))= |G|      : True
        All i <> j,            Sum[g](Xi(g)Xj(g^−1))=  0       : True
        All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|  : True
        All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0       : True

        |Symm6| = 720
        Type        NonAbelianGroup
        BaseGroup   S6

        [Class      1 2a 2b 2c 3a 3b 4a 4b   5  6a  6b]
        [ Size      1 15 15 45 40 40 90 90 144 120 120]
        [                                             ]
        [  X.1      1  1  1  1  1  1  1  1   1   1   1]
        [  X.2      1 -1 -1  1  1  1 -1  1   1  -1  -1]
        [  X.3      5  3 -1  1  2 -1  1 -1   0   0  -1]
        [  X.4      5 -3  1  1  2 -1 -1 -1   0   0   1]
        [  X.5      5 -1  3  1 -1  2  1 -1   0  -1   0]
        [  X.6      5  1 -3  1 -1  2 -1 -1   0   1   0]
        [  X.7      9  3  3  1  0  0 -1  1  -1   0   0]
        [  X.8      9 -3 -3  1  0  0  1  1  -1   0   0]
        [  X.9     10  2 -2 -2  1  1  0  0   0  -1   1]
        [ X.10     10 -2  2 -2  1  1  0  0   0   1  -1]
        [ X.11     16  0  0  0 -2 -2  0  0   1   0   0]
        All i,                 Sum[g](Xi(g)Xi(g^−1))= |G|      : True
        All i <> j,            Sum[g](Xi(g)Xj(g^−1))=  0       : True
        All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|  : True
        All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0       : True
     */

    public static void ExamplesGL23()
    {
        var gl = new GL(2, 3);
        var a0 = gl[2, 0, 0, 1];
        var b0 = gl[2, 1, 2, 0];

        var gl23 = Group.Generate(gl, a0, b0);
        var ctGL23 = FG.CharactersTable2Slow(gl23);
        ctGL23.DisplayCells();
    
        var a = gl[1, 1, 0, 1];
        var b = gl[0, 1, 2, 0];

        var sl23 = Group.Generate("SL(2,3)", gl, a, b);
        var ctSL23 = FG.CharactersTable2(sl23);
        ctSL23.RestrictionFromSuperGroup(ctGL23);
        ctSL23.DisplayCells();
    }
    /*
        |GL(2,3)| = 48
        Type        NonAbelianGroup
        BaseGroup   GL(2,3)

        [Class     1 2a 2b  3  4  6         8a         8b]
        [ Size     1  1 12  8  6  8          6          6]
        [                                                ]
        [  X.1     1  1  1  1  1  1          1          1]
        [  X.2     1  1 -1  1  1  1         -1         -1]
        [  X.3     2  2  0 -1  2 -1          0          0]
        [  X.4     2 -2  0 -1  0  1   ξ8³ + ξ8 -ξ8³ + -ξ8]
        [  X.5     2 -2  0 -1  0  1 -ξ8³ + -ξ8   ξ8³ + ξ8]
        [  X.6     3  3 -1  0 -1  0          1          1]
        [  X.7     3  3  1  0 -1  0         -1         -1]
        [  X.8     4 -4  0  1  0 -1          0          0]
        All i,                 Sum[g](Xi(g)Xi(g^−1))= |G|      : True
        All i <> j,            Sum[g](Xi(g)Xj(g^−1))=  0       : True
        All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|  : True
        All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0       : True

        |SL(2,3)| = 24
        Type        NonAbelianGroup
        BaseGroup   GL(2,3)

        [Class     1  2       3a       3b  4       6a       6b]
        [ Size     1  1        4        4  6        4        4]
        [                                                     ]
        [  X.1     1  1        1        1  1        1        1]
        [  X.2     1  1 -ξ3 + -1       ξ3  1       ξ3 -ξ3 + -1]
        [  X.3     1  1       ξ3 -ξ3 + -1  1 -ξ3 + -1       ξ3]
        [  X.4     2 -2       -1       -1  0        1        1]
        [  X.5     2 -2   ξ3 + 1      -ξ3  0       ξ3 -ξ3 + -1]
        [  X.6     2 -2      -ξ3   ξ3 + 1  0 -ξ3 + -1       ξ3]
        [  X.7     3  3        0        0 -1        0        0]
        All i,                 Sum[g](Xi(g)Xi(g^−1))= |G|      : True
        All i <> j,            Sum[g](Xi(g)Xj(g^−1))=  0       : True
        All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1))= |Cl(g)|  : True
        All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1))=  0       : True
     */
    
    public static void OtherExamples()
    {
        {
            FG.CharactersTable2(FG.DiCyclicSdp(5)).DisplayCells();
            FG.CharactersTable2(FG.DiCyclicSdp(6)).DisplayCells();
            FG.CharactersTable2(FG.DiCyclicSdp(7)).DisplayCells();
            FG.CharactersTable2(FG.DiCyclicSdp(8)).DisplayCells();
        }

        {
            FG.CharactersTable2(Group.SemiDirectProd(FG.Abelian(4, 4), FG.Abelian(3))).DisplayCells();
            FG.CharactersTable2(Group.SemiDirectProd(FG.Abelian(3), FG.Symmetric(3))).DisplayCells();
            FG.CharactersTable2(Group.SemiDirectProd(FG.Abelian(3, 3), FG.Abelian(8))).DisplayCells();
            FG.CharactersTable2(Group.SemiDirectProd(FG.Abelian(5), FG.Abelian(8))).DisplayCells();
            FG.CharactersTable2(Group.SemiDirectProd(FG.Abelian(4), FG.Abelian(4))).DisplayCells();
            FG.CharactersTable2(Group.SemiDirectProd(FG.Abelian(4), FG.Abelian(4))).DisplayCells();
        }

        {
            var E = FG.WordGroup("E", "a4, d6, adad, a3da3d");
            var Esdp = Group.SemiDirectProd(FG.Abelian(3), FG.DihedralSdp(4));
            DisplayGroup.AreIsomorphics(Esdp, E);
            FG.CharactersTable2Slow(E).DisplayCells();
            FG.CharactersTable2Slow(Esdp).DisplayCells();
        }
        
        {
            var f8 = Group.SemiDirectProd(FG.Abelian(2, 2, 2), FG.Abelian(7));
            FG.CharactersTable2(f8).DisplayCells();
            
            var gr = Group.SemiDirectProd(f8, FG.Abelian(3));
            FG.CharactersTable2Slow(gr).DisplayCells();
        }
    }
}