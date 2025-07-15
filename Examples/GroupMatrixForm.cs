using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Matrix;

namespace Examples;

public static class GroupMatrixForm
{
    static (int m, int n, int r)[] MetaCyclicSdp(int order)
    {
        return IntExt.Dividors(order).Where(d => d > 1)
            .SelectMany(m => FG.MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
            .ToArray();
    }

    static ConcreteGroup<Mat> MetaCyclicGL2p_Meth2(int m, int n, int r)
    {
        foreach (var p in IntExt.Primes10000.Where(p => m % p == 0 && (p - 1) % n == 0))
        {
            var Fp = FG.UnInt(p);
            var gl = new GL(2, p);

            var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0).ToArray();
            var ordn = Fp.Where(e => Fp.ElementsOrders[e] == n).ToArray();

            var m0s = ordm.Select(e => gl[e.K, 1, 0, e.K])
                .Where(mat => mat.IsOrder(m))
                .ToArray();

            var m1s = ordn.Grid2D(ordn.Append(Fp.Neutral())).Select(e => gl[e.t1.K, 0, 0, e.t2.K])
                .Where(mat => mat.IsOrder(n))
                .ToArray();

            foreach (var (m0, m1) in m0s.Grid2D(m1s).Where(e => gl.Op(gl.Invert(e.t2), gl.Op(e.t1, e.t2)).Equals(gl.Times(e.t1, r))))
            {
                var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
                if (mtGl.Count() == m * n)
                    return mtGl;
            }
        }

        return Group.Generate(new GL(1, 2));
    }

    static void AllGensOfMtCycSdpUpToOrder(int maxOrd, int maxDim = 6, bool altGL2Meth = true)
    {
        GlobalStopWatch.Restart();
        var missing = new List<(int, int, int)>();
        var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => MetaCyclicSdp(ord)).ToArray();

        foreach (var e in allMtCycSdp)
        {
            if (altGL2Meth)
            {
                var mtGL2 = MetaCyclicGL2p_Meth2(e.m, e.n, e.r);
                if (mtGL2.Count() != 1)
                {
                    DisplayGroup.HeadOrdersGenerators(mtGL2);
                    continue;
                }
            }
            
            var mtGL = FG.MetaCyclicSdpMat(e.m, e.n, e.r, maxDim);
            if (mtGL.Count() != 1)
            {
                DisplayGroup.HeadOrdersGenerators(mtGL);
                continue;
            }

            missing.Add(e);
        }

        var total = allMtCycSdp.Length;
        missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
        GlobalStopWatch.Show("END");
        Console.Beep();
    }
    
    public static void ExampleDihedalGL2p()
    {
        for (int n = 3; n < 33; n++)
        {
            var D2n = FG.DihedralGL2p(n);
            var D2pg = FG.Dihedral(n);
            D2pg.Name = $"{D2pg}pg";

            DisplayGroup.Generators(D2n);
            DisplayGroup.Generators(D2pg);

            if (!D2n.IsIsomorphicTo(D2pg))
                throw new();

            Console.WriteLine($"{D2n} IsIsomorphicTo {D2pg}");
            Console.WriteLine();
        }
    }
    
    public static void ExampleAbelian()
    {
        GlobalStopWatch.Restart();
        var total = 0;
        var maxP = 0;
        for (int n = 1; n <= 256; n++)
        {
            var allAb = FG.AllAbelianGroupsOfOrder(n);
            Console.WriteLine($"############ Abelian groups or Order {n} ############");
            foreach (var ab in allAb)
            {
                ++total;
                var seq = ab.Neutral().Ei.Select(e => e.Mod).ToArray();
                var abMat = FG.AbelianMat(seq);

                DisplayGroup.HeadGenerators(abMat);
                if (!ab.IsIsomorphicTo(abMat))
                    throw new();

                maxP = int.Max(maxP, abMat.Neutral().GL.P);
                Console.WriteLine();
            }
        }

        GlobalStopWatch.Show($"Total:{total} Max Prime:{maxP}");
    }
    
    public static void ExampleDicyclic()
    {
        for (int m = 2; m < 33; m++)
        {
            var Dic_m = FG.DicyclicGL2p(m);
            var Dic_m_wg = FG.DiCyclic(m);
            Dic_m_wg.Name = $"{Dic_m_wg}_wg";
            DisplayGroup.HeadGenerators(Dic_m);
            if (!Dic_m.IsIsomorphicTo(Dic_m_wg))
                throw new();

            Console.WriteLine($"{Dic_m} IsIsomorphicTo {Dic_m_wg}");
            Console.WriteLine();
        }
    }
    
    public static void ExampleSemiDihedralAndModularMax()
    {
        for (int n = 3; n < 9; n++)
        {
            var qd = FG.SemiDihedralGL2p(n);
            DisplayGroup.HeadGenerators(qd);
            if (!qd.IsIsomorphicTo(FG.SemiDihedral(n)))
                throw new();

            var mm = FG.ModularMaxGL2p(n);
            DisplayGroup.HeadGenerators(mm);
            if (!mm.IsIsomorphicTo(FG.ModularMax(n)))
                throw new();

            Console.WriteLine();
        }
    }

    public static void ExampleAllMetaCyclicSemiDirectProducts()
    {
        // Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracketNoFmt;
        Group.ActivedStorage(false);

        AllGensOfMtCycSdpUpToOrder(32);
        // AllGensOfMtCycSdpUpToOrder(64);
        // AllGensOfMtCycSdpUpToOrder(maxOrd: 256, maxDim: 12, altGL2Meth: false);
    }
}