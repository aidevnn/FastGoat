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

    public static void ExamplesCharacterTableAbelianGroups()
    {
        FG.CharacterTable(FG.Abelian(3)).DisplayCells();
        FG.CharacterTable(FG.Abelian(2, 2)).DisplayCells();
        FG.CharacterTable(FG.Abelian(2, 3)).DisplayCells();
        FG.CharacterTable(FG.Abelian(6)).DisplayCells();
        FG.CharacterTable(FG.Abelian(4)).DisplayCells();
        FG.CharacterTable(FG.Abelian(2, 4)).DisplayCells();
        FG.CharacterTable(FG.Abelian(2, 2, 2)).DisplayCells();
        FG.CharacterTable(FG.Abelian(2, 2, 3)).DisplayCells();
    }

    public static void MoreExamples()
    {
        FG.CharacterTable(FG.Symmetric(3)).DisplayCells();
        FG.CharacterTable(FG.Dihedral(4)).DisplayCells();
        FG.CharacterTable(FG.Alternate(4)).DisplayCells();
        FG.CharacterTable(FG.Quaternion(8)).DisplayCells();
        FG.CharacterTable(FG.Dihedral(5)).DisplayCells();
        FG.CharacterTable(FG.Dihedral(8)).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(5), new Cn(4))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(7), new Cn(6))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(11), new Cn(10))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(13), new Cn(12))).DisplayCells();

        FG.CharacterTable(Group.SemiDirectProd(new Cn(3), new Cn(4))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(4), new Cn(4))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(7), new Cn(3))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(9), new Cn(3))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(new Cn(3), new Cn(8))).DisplayCells();
        FG.CharacterTable(Group.SemiDirectProd(FG.Abelian(3, 3), new Cn(4))).DisplayCells();

        FG.CharacterTable(FG.Dihedral(6)).DisplayCells();
        FG.CharacterTable(FG.DiCyclic(3)).DisplayCells();
        FG.CharacterTable(FG.SemiDihedral(4)).DisplayCells();
    }

    public static void ExamplesPQgroups()
    {
        for (int i = 2; i <= 16; i++)
            FG.CharacterTable(FG.DihedralSdp(i)).DisplayCells();

        var pqGroups = 31.Range(2).SelectMany(i => FG.FrobeniusSdp(i)).ToArray();
        foreach (var g in pqGroups)
        {
            FG.CharacterTable(g).DisplayCells();
        }
    }

    public static void ExamplesPermGroups()
    {
        GlobalStopWatch.Restart();
        for (int n = 3; n < 7; n++)
        {
            GlobalStopWatch.AddLap();
            var ctAn = FG.CharacterTable(FG.Alternate(n));
            ctAn.DisplayCells();
            GlobalStopWatch.Show($"{ctAn.Gr}");
        }

        for (int n = 3; n < 7; n++)
        {
            GlobalStopWatch.AddLap();
            var ctSn = FG.CharacterTable(FG.Symmetric(n));
            ctSn.DisplayCells();
            GlobalStopWatch.Show($"{ctSn.Gr}");
        }
    }
    /*
       |Alt5| = 60
       Type        NonAbelianGroup
       BaseGroup   S5
       SuperGroup  |Symm5| = 120

       [Class      1   2   3            5a            5b]
       [ Size      1  15  20            12            12]
       [                                                ]
       [  Ꭓ.1      1   1   1             1             1]
       [  Ꭓ.2      3  -1   0  1/2 - 1/2·√5  1/2 + 1/2·√5]
       [  Ꭓ.3      3  -1   0  1/2 + 1/2·√5  1/2 - 1/2·√5]
       [  Ꭓ.4      4   0   1            -1            -1]
       [  Ꭓ.5      5   1  -1             0             0]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True

       |Symm5| = 120
       Type        NonAbelianGroup
       BaseGroup   S5

       [Class      1  2a  2b   3   4   5   6]
       [ Size      1  10  15  20  30  24  20]
       [                                    ]
       [  Ꭓ.1      1   1   1   1   1   1   1]
       [  Ꭓ.2      1  -1   1   1  -1   1  -1]
       [  Ꭓ.3      4   2   0   1   0  -1  -1]
       [  Ꭓ.4      4  -2   0   1   0  -1   1]
       [  Ꭓ.5      5   1   1  -1  -1   0   1]
       [  Ꭓ.6      5  -1   1  -1   1   0  -1]
       [  Ꭓ.7      6   0  -2   0   0   1   0]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True

       |Alt6| = 360
       Type        NonAbelianGroup
       BaseGroup   S6
       SuperGroup  |Symm6| = 720

       [Class       1   2  3a  3b   4            5a            5b]
       [ Size       1  45  40  40  90            72            72]
       [                                                         ]
       [  Ꭓ.1       1   1   1   1   1             1             1]
       [  Ꭓ.2       5   1  -1   2  -1             0             0]
       [  Ꭓ.3       5   1   2  -1  -1             0             0]
       [  Ꭓ.4       8   0  -1  -1   0  1/2 - 1/2·√5  1/2 + 1/2·√5]
       [  Ꭓ.5       8   0  -1  -1   0  1/2 + 1/2·√5  1/2 - 1/2·√5]
       [  Ꭓ.6       9   1   0   0   1            -1            -1]
       [  Ꭓ.7      10  -2   1   1   0             0             0]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True

       |Symm6| = 720
       Type        NonAbelianGroup
       BaseGroup   S6

       [Class       1  2a  2b  2c  3a  3b  4a  4b   5  6a  6b]
       [ Size       1  15  15  45  40  40  90  90 144 120 120]
       [                                                     ]
       [  Ꭓ.1       1   1   1   1   1   1   1   1   1   1   1]
       [  Ꭓ.2       1  -1  -1   1   1   1  -1   1   1  -1  -1]
       [  Ꭓ.3       5   1  -3   1  -1   2  -1  -1   0   1   0]
       [  Ꭓ.4       5  -1   3   1  -1   2   1  -1   0  -1   0]
       [  Ꭓ.5       5   3  -1   1   2  -1   1  -1   0   0  -1]
       [  Ꭓ.6       5  -3   1   1   2  -1  -1  -1   0   0   1]
       [  Ꭓ.7       9   3   3   1   0   0  -1   1  -1   0   0]
       [  Ꭓ.8       9  -3  -3   1   0   0   1   1  -1   0   0]
       [  Ꭓ.9      10   2  -2  -2   1   1   0   0   0  -1   1]
       [ Ꭓ.10      10  -2   2  -2   1   1   0   0   0   1  -1]
       [ Ꭓ.11      16   0   0   0  -2  -2   0   0   1   0   0]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True
     */

    public static void ExamplesGL23()
    {
        var gl = new GL(2, 3);
        var a0 = gl[2, 0, 0, 1];
        var b0 = gl[2, 1, 2, 0];

        var gl23 = Group.Generate(gl, a0, b0);
        var ctGL23 = FG.CharacterTable(gl23);
        ctGL23.DisplayCells();

        var a = gl[1, 1, 0, 1];
        var b = gl[0, 1, 2, 0];

        var sl23 = Group.Generate("SL(2,3)", gl, a, b);
        var ctSL23 = FG.CharacterTableEmpty(sl23);
        ctSL23.DerivedSubGroupLift();
        ctSL23.DisplayCells();
        ctSL23.RestrictionFromSuperGroup(ctGL23);
        ctSL23.DisplayCells();
    }
    /*
       |GL(2,3)| = 48
       Type        NonAbelianGroup
       BaseGroup   GL(2,3)

       [Class      1  2a  2b   3   4   6    8a    8b]
       [ Size      1   1  12   8   6   8     6     6]
       [                                            ]
       [  Ꭓ.1      1   1   1   1   1   1     1     1]
       [  Ꭓ.2      1   1  -1   1   1   1    -1    -1]
       [  Ꭓ.3      2   2   0  -1   2  -1     0     0]
       [  Ꭓ.4      2  -2   0  -1   0   1  -I√2   I√2]
       [  Ꭓ.5      2  -2   0  -1   0   1   I√2  -I√2]
       [  Ꭓ.6      3   3   1   0  -1   0    -1    -1]
       [  Ꭓ.7      3   3  -1   0  -1   0     1     1]
       [  Ꭓ.8      4  -4   0   1   0  -1     0     0]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True

       |SL(2,3)| = 24
       Type        NonAbelianGroup
       BaseGroup   GL(2,3)

       [Class      1   2    3a    3b   4   6a   6b]
       [ Size      1   1     4     4   6    4    4]
       [                                          ]
       [  Ꭓ.1      1   1     1     1   1    1    1]
       [  Ꭓ.2      1   1    ξ3   ξ3²   1  ξ3²   ξ3]
       [  Ꭓ.3      1   1   ξ3²    ξ3   1   ξ3  ξ3²]
       [  Ꭓ.4      2  -2    -1    -1   0    1    1]
       [  Ꭓ.5      2  -2   -ξ3  -ξ3²   0  ξ3²   ξ3]
       [  Ꭓ.6      2  -2  -ξ3²   -ξ3   0   ξ3  ξ3²]
       [  Ꭓ.7      3   3     0     0  -1    0    0]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True
     */

    public static void OtherExamples()
    {
        {
            FG.CharacterTable(FG.DiCyclicSdp(5)).DisplayCells();
            FG.CharacterTable(FG.DiCyclicSdp(6)).DisplayCells();
            FG.CharacterTable(FG.DiCyclicSdp(7)).DisplayCells();
            FG.CharacterTable(FG.DiCyclicSdp(8)).DisplayCells();
        }

        {
            FG.CharacterTable(Group.SemiDirectProd(FG.Abelian(4, 4), FG.Abelian(3))).DisplayCells();
            FG.CharacterTable(Group.SemiDirectProd(FG.Abelian(3), FG.Symmetric(3))).DisplayCells();
            FG.CharacterTable(Group.SemiDirectProd(FG.Abelian(3, 3), FG.Abelian(8))).DisplayCells();
            FG.CharacterTable(Group.SemiDirectProd(FG.Abelian(5), FG.Abelian(8))).DisplayCells();
            FG.CharacterTable(Group.SemiDirectProd(FG.Abelian(4), FG.Abelian(4))).DisplayCells();
        }

        {
            var E = FG.WordGroup("E", "a4, d6, adad, a3da3d");
            var Esdp = Group.SemiDirectProd(FG.Abelian(3), FG.DihedralSdp(4));
            DisplayGroup.AreIsomorphics(Esdp, E);
            FG.CharacterTable(E).DisplayCells();
            FG.CharacterTable(Esdp).DisplayCells();
        }

        {
            var f8 = Group.SemiDirectProd(FG.Abelian(2, 2, 2), FG.Abelian(7));
            FG.CharacterTable(f8).DisplayCells();

            var gr = Group.SemiDirectProd(f8, FG.Abelian(3));
            FG.CharacterTable(gr).DisplayCells();
        }
    }

    public static void Symmetric7CharacterTable()
    {
        GlobalStopWatch.Restart();

        GlobalStopWatch.AddLap();
        var ctA7 = FG.CharacterTableEmpty(FG.Alternate(7));
        ctA7.DisplayCells();

        ctA7.InductionFromSubGroup(Group.Generate("Alt6", ctA7.Gr, ctA7.Gr[(1, 2, 3)], ctA7.Gr[(1, 2, 4)], ctA7.Gr[(1, 2, 5)],
            ctA7.Gr[(1, 2, 6)]));
        ctA7.DisplayCells();

        ctA7.InductionFromSubGroup(Group.Generate("S5", ctA7.Gr, ctA7.Gr[(1, 2), (3, 4), (5, 6, 7)], ctA7.Gr[(3, 4, 5, 6, 7)]));
        ctA7.DisplayCells();

        ctA7.InductionFromSubGroup(Group.Generate("C3 x A4", ctA7.Gr, ctA7.Gr[(1, 2), (3, 4), (5, 6, 7)],
            ctA7.Gr[(2, 3, 4), (5, 6, 7)]));
        ctA7.DisplayCells();

        ctA7.InductionFromSubGroup(Group.Generate("L(7)", ctA7.Gr, ctA7.Gr[(1, 2, 3, 4, 5, 6, 7)], ctA7.Gr[(1, 2), (3, 6)]));
        ctA7.DisplayCells();
        GlobalStopWatch.Show($"{ctA7.Gr}");

        GlobalStopWatch.AddLap();
        var ctS7 = FG.CharacterTableEmpty(FG.Symmetric(7));
        ctS7.DisplayCells();

        ctS7.InductionFromSubGroup(ctA7);
        ctS7.DisplayCells();

        ctS7.InductionFromSubGroup(Group.Generate("Symm6", ctS7.Gr, ctS7.Gr[(1, 2)], ctS7.Gr[(1, 2, 3, 4, 5, 6)]));
        ctS7.DisplayCells();

        ctS7.InductionFromSubGroup(Group.Generate("C2 x S4", ctS7.Gr, ctS7.Gr[(6, 7)], ctS7.Gr[(4, 5)], ctS7.Gr[(3, 4)],
            ctS7.Gr[(2, 3)],
            ctS7.Gr[(1, 2)]));
        ctS7.DisplayCells();
        GlobalStopWatch.Show($"{ctS7.Gr}");
    }

    /*
       |Alt7| = 2520
       Type        NonAbelianGroup
       BaseGroup   S7
       SuperGroup  |Symm7| = 5040

       [Class       1   2  3a  3b   4   5   6              7a              7b]
       [ Size       1 105  70 280 630 504 210             360             360]
       [                                                                     ]
       [  Ꭓ.1       1   1   1   1   1   1   1               1               1]
       [  Ꭓ.2       6   2   3   0   0   1  -1              -1              -1]
       [  Ꭓ.3      10  -2   1   1   0   0   1  -1/2 - 1/2·I√7  -1/2 + 1/2·I√7]
       [  Ꭓ.4      10  -2   1   1   0   0   1  -1/2 + 1/2·I√7  -1/2 - 1/2·I√7]
       [  Ꭓ.5      14   2  -1   2   0  -1  -1               0               0]
       [  Ꭓ.6      14   2   2  -1   0  -1   2               0               0]
       [  Ꭓ.7      15  -1   3   0  -1   0  -1               1               1]
       [  Ꭓ.8      21   1  -3   0  -1   1   1               0               0]
       [  Ꭓ.9      35  -1  -1  -1   1   0  -1               0               0]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True

       # Alt7 Time:92077 ms

       |Symm7| = 5040
       Type        NonAbelianGroup
       BaseGroup   S7

       [Class       1  2a  2b  2c  3a  3b  4a  4b   5  6a  6b  6c   7  10  12]
       [ Size       1  21 105 105  70 280 210 630 504 210 420 840 720 504 420]
       [                                                                     ]
       [  Ꭓ.1       1   1   1   1   1   1   1   1   1   1   1   1   1   1   1]
       [  Ꭓ.2       1  -1   1  -1   1   1  -1   1   1   1  -1  -1   1  -1  -1]
       [  Ꭓ.3       6   4   2   0   3   0   2   0   1  -1   1   0  -1  -1  -1]
       [  Ꭓ.4       6  -4   2   0   3   0  -2   0   1  -1  -1   0  -1   1   1]
       [  Ꭓ.5      14   4   2   0  -1   2  -2   0  -1  -1   1   0   0  -1   1]
       [  Ꭓ.6      14  -4   2   0  -1   2   2   0  -1  -1  -1   0   0   1  -1]
       [  Ꭓ.7      14   6   2   2   2  -1   0   0  -1   2   0  -1   0   1   0]
       [  Ꭓ.8      14  -6   2  -2   2  -1   0   0  -1   2   0   1   0  -1   0]
       [  Ꭓ.9      15   5  -1  -3   3   0   1  -1   0  -1  -1   0   1   0   1]
       [ Ꭓ.10      15  -5  -1   3   3   0  -1  -1   0  -1   1   0   1   0  -1]
       [ Ꭓ.11      20   0  -4   0   2   2   0   0   0   2   0   0  -1   0   0]
       [ Ꭓ.12      21   1   1  -3  -3   0  -1  -1   1   1   1   0   0   1  -1]
       [ Ꭓ.13      21  -1   1   3  -3   0   1  -1   1   1  -1   0   0  -1   1]
       [ Ꭓ.14      35   5  -1   1  -1  -1  -1   1   0  -1  -1   1   0   0  -1]
       [ Ꭓ.15      35  -5  -1  -1  -1  -1   1   1   0  -1   1  -1   0   0   1]
       All i,                 Sum[g](Xi(g)Xi(g^−1)) = |G|      : True
       All i <> j,            Sum[g](Xi(g)Xj(g^−1)) =  0       : True
       All g, h in Cl(g),     Sum[r](Xr(g)Xr(h^−1)) = |Cl(g)|  : True
       All g, h not in Cl(g), Sum[r](Xr(g)Xr(h^−1)) =  0       : True

       # Symm7 Time:379028 ms
     */
}