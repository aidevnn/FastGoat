using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;

namespace FastGoat.Examples;

public static class GroupMatrixForm
{
    static ConcreteGroup<Mat> AbelianGroup(params int[] seq)
    {
        var dim = seq.Length;
        if (dim > 5)
            throw new("Dim 5 Maximum Matrix");
        
        var p = IntExt.Primes10000.First(p => seq.All(o => (p - 1) % o == 0));
        var gl = new GL(dim, p);
        var a0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
        var a = new ZnInt(p, a0);
        var Zp = Group.MulGroup($"Z({p})", a);
        var seq2 = seq.Select(o => Zp.ElementsOrders.First(e => e.Value == o).Key.K).ToArray();

        int[] Diag(int n, int k, int v) => n.Range().Grid2D()
            .Select(e => e.t1 == e.t2 ? (e.t1 == k ? v : 1) : 0).ToArray();

        var gens = seq2.Select((v, k) => gl.Create(Diag(dim, k, v))).ToArray();
        return Group.Generate(seq.Glue(" x ", "C{0}"), gl, gens);
    }

    static int OrderMatOrth(ZnInt x0, ZnInt y0)
    {
        var (x1, y1) = (x0.One, y0.Zero);
        for (int i = 1; i < 1000; i++)
        {
            (x1, y1) = (x0 * x1 - y0 * y1, x0 * y1 + y0 * x1);
            if (x1.Equals(x1.One) && y1.IsZero())
                return i;
        }

        throw new("####################################");
    }

    static ConcreteGroup<Mat> DihedralGO2p(int n = 16)
    {
        int p = 2;
        ZnInt a, x = ZnInt.ZnZero(p), y = ZnInt.ZnZero(p);
        foreach (var p0 in IntExt.Primes10000.Where(p0 => (p0 - 1) % n == 0))
        {
            p = p0;
            var a0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
            a = new ZnInt(p, a0);
            var Zp = Group.MulGroup($"F{p}", a);
            var square = Zp.Append(a.Zero).Select(x0 => (x:x0, x2: x0 * x0)).GroupBy(e => e.x2)
                .ToDictionary(e => e.Key, e => e.Select(f => f.x).ToArray());
            var dicSquare = Zp.ToDictionary(x0 => x0, x0 => square.ContainsKey(x0) ? square[x0] : []);
            dicSquare[a.Zero] = [];
            var XYs = Zp.Append(a.Zero)
                .Select(x0 => (x: x0, yList: dicSquare[1 - x0 * x0]))
                .Where(e => e.yList.Length != 0)
                .SelectMany(e => e.yList.Select(y0 => (e.x, y: y0)))
                .Distinct()
                .Select(e => (e.x, e.y))
                .Where(e => OrderMatOrth(e.x, e.y) == n)
                .OrderBy(e => e.x.K)
                .ToArray();

            if (XYs.Length != 0)
            {
                (x, y) = XYs[0];
                break;
            }
        }

        var gl = new GL(2, p);
        var m0 = gl[x.K, y.K, (-y).K, x.K];
        var m1 = gl[0, 1, 1, 0];

        return Group.Generate($"D{2 * n}", gl, m0, m1);
    }

    public static void ExampleDihedalGO2p()
    {
        for (int n = 3; n < 33; n++)
        {
            var D2n = DihedralGO2p(n);
            var D2pg = FG.Dihedral(n);
            D2pg.Name = $"{D2pg}pg";
            
            DisplayGroup.Generators(D2n);
            DisplayGroup.Generators(D2pg);
            
            DisplayGroup.AreIsomorphics(D2n, D2pg);
            Console.WriteLine();
        }
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
        var maxP = 1;
        for (int n = 3; n < 64; n++)
        {
            var allAb = FG.AllAbelianGroupsOfOrder(n);
            Console.WriteLine($"############ Abelian groups or Order {n} ############");
            foreach (var ab in allAb)
            {
                var seq = ab.Neutral().Ei.Select(e => e.Mod).ToArray();
                var abMat = AbelianGroup(seq);
                
                DisplayGroup.Generators(abMat);
                DisplayGroup.Generators(ab);
                if (!ab.IsIsomorphicTo(abMat))
                    throw new();
            
                Console.WriteLine($"{abMat} IsIsomorphicTo {ab}");
                Console.WriteLine();

                maxP = int.Max(maxP, abMat.Neutral().GL.P);
            }
            
            Console.WriteLine();
        }

        Console.WriteLine(new { maxP });
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

}