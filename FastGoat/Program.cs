using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions.Group;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

bool IsOrder(KMatrix<ZnInt> m, int o)
{
    var e = m.One;
    var mk = m.One;
    if (!m.Det.Pow(o).Equals(m.KOne))
        return false;

    for (int k = 0; k < o; k++)
    {
        if (k != 0 && mk.Equals(e))
            return false;

        mk *= m;
    }

    return mk.Equals(e);
}

int[][] DistinctPartition(int dim) => Partitions32[dim].Where(l => l.Count == l.Distinct().Count())
    .Select(l => l.Order().ToArray())
    .OrderBy(l => l.Length)
    .ToArray();

IEnumerable<(int[] perm, int[][] cycles)> PartitionsToPerm(int dim)
{
    var list = new List<(int[] perm, int[][] cycles)>();
    foreach (var type in DistinctPartition(dim))
    {
        var lt = new List<int[]>();
        var rg = dim.Range();
        var perm = new int[dim];
        foreach (var k in type)
        {
            var r0 = rg.Take(k).ToArray();
            rg = rg.Skip(k).ToArray();
            lt.Add(r0);
            for (int i = 0; i < k; i++)
                perm[r0[i]] = r0[(i + 1) % k];
        }

        yield return (perm, lt.ToArray());
    }
}

IEnumerable<(int[] perm, int[][] cycles, KMatrix<ZnInt> mat)> PartitionsToMatrixOfOrder(GLn<ZnInt> GL, int ord)
{
    var dim = GL.N;
    var sn = new Sn(dim);
    var o = GL.Neutral().Rows.Select(r => r.ToArray()).ToArray();
    
    var Fp = FG.UnInt(GL.Neutral().P);
    var ordn = Fp.Where(e => ord % Fp.ElementsOrders[e] == 0)
        .Select(e =>
        {
            var diag = Ring.Diagonal(e.One, dim);
            diag[0, 0] = e;
            return new KMatrix<ZnInt>(diag);
        }).ToArray();
    
    foreach (var (perm, cycles) in PartitionsToPerm(dim))
    {
        var perm0 = sn.CreateElement(perm.Select(i => i + 1).ToArray());
        var o0 = perm0.Apply(o.Select(l => l.ToArray()).ToArray());
        var mat0 = o0.SelectMany(r => r).ToKMatrix(dim);
        foreach (var m in ordn)
        {
            var mat = m * mat0;
            if (IsOrder(mat, ord))
                yield return (perm, cycles, mat);
        }
    }
}

ConcreteGroup<KMatrix<ZnInt>> MetaCyclicGLnp_DiagByPerm(int m, int n, int r, int dim)
{
    var nks = DistinctPartition(dim).Select(l => l.Aggregate((a0, a1) => a0 * a1))
        .SelectMany(e => Dividors(e).Append(e).Where(j => j != 1))
        .Append(n)
        .Distinct()
        .ToArray();
    foreach (var p in nks.Select(nk => Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % nk == 0)).Distinct().Order())
    {
        var Fp = FG.UnInt(p);
        var o = Fp.Neutral();
        var GL = FG.GLnK($"F{p}", dim, o);
        // Console.WriteLine($"Solve M({m}x:{n}){r} in {GL}");

        foreach (var (perm, cycles, m1) in PartitionsToMatrixOfOrder(GL, n))
        {
            var m1i = m1.Inv();
            var seq = cycles.Select(c => c.Length).Select(l =>
            {
                var r0 = IntExt.PowMod(r, l, m);
                var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0 && e.Pow(r0).Equals(e))
                    .OrderByDescending(e => Fp.ElementsOrders[e]).ToArray();

                return ordm.Select(a => l.Range(1).Select(k => a.Pow(IntExt.PowMod(r, k, m))).Reverse().ToArray());
            }).MultiLoop().Select(l => l.ToArray()).Where(l => l.All(l1 => l1.Any(e => !e.Equals(e.One))));
            foreach (var l in seq.Take(1))
            {
                var arr = new int[dim * dim];
                foreach (var (sols, idxs) in l.Zip(cycles))
                foreach (var (idx, sol) in idxs.Zip(sols))
                    arr[perm[idx] * (dim + 1)] = sol.K;

                var m0 = arr.Select(i => i * o).ToKMatrix(dim);
                if (IsOrder(m0, m) && (m1i * m0 * m1).Equals(m0.Pow(r)))
                {
                    var mtGL = Group.Generate($"M({m}x:{n}){r}", GL, m0, m1);
                    if (mtGL.Count() == m * n)
                    {
                        Console.WriteLine($"{mtGL} Perm type[{cycles.Select(l0 => l0.Length).Glue(",")}] rank:{cycles.Length}");
                        return mtGL;
                    }
                }
            }
        }
    }

    return Group.Generate(FG.GLnK("F2", 1, ZnInt.ZnZero(2)));
}

(int m, int n, int r)[] MetaCyclicSdp(int order)
{
    return IntExt.Dividors(order).Where(d => d > 1)
        .SelectMany(m => FG.MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
        .ToArray();
}

void AllGensOfMtCycSdpUpToOrder(int maxOrd, int maxDim)
{
    GlobalStopWatch.Restart();
    var missing = new List<(int, int, int)>();
    var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => MetaCyclicSdp(ord)).ToArray();

    foreach (var e in allMtCycSdp)
    {
        var found = false;
        foreach (var dim in maxDim.Range(1).Where(d => d != 1 && (Gcd(e.m, d) != 1 || Gcd(e.m - 1, d) != 1)))
        {
            var mtGL = MetaCyclicGLnp_DiagByPerm(e.m, e.n, e.r, dim);
            if (mtGL.Count() != 1)
            {
                found = true;
                DisplayGroup.HeadOrdersGenerators(mtGL);
                break;
            }
        }

        if (!found)
            missing.Add(e);
    }

    var total = allMtCycSdp.Length;
    missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    GlobalStopWatch.Show("End Gens");
    Console.Beep();
}

{
    // Ring.MatrixDisplayForm = Ring.MatrixDisplay.OneLineArray;
    
    AllGensOfMtCycSdpUpToOrder(maxOrd: 256, maxDim: 12);
}