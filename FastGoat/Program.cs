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

(int[] perm, int[][] cycles)[] PartitionsToPerm(int dim)
{
    var all = Partitions32[dim].Select(l => l.Order().ToArray()).ToArray();
    var list = new List<(int[] perm, int[][] cycles)>();
    foreach (var type in all)
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

        list.Add((perm, lt.ToArray()));
    }

    return list.ToArray();
}

(int[] perm, int[][] cycles, KMatrix<ZnInt> mat)[] PartitionsToMatrix(GLn<ZnInt> GL, ZnInt e)
{
    var dim = GL.N;
    var sn = new Sn(dim);
    var seq = PartitionsToPerm(dim);
    var o = GL.Neutral().Rows.Select(r => r.ToArray()).ToArray();
    var list = new List<(int[] perm, int[][] cycles, KMatrix<ZnInt>)>();
    foreach (var (perm, cycles) in seq)
    {
        var perm0 = sn.CreateElement(perm.Select(i => i + 1).ToArray());
        var o0 = perm0.Apply(o.Select(l => l.ToArray()).ToArray());
        o0[0] = o0[0].Select(c => c * e).ToArray();
        var mat = o0.SelectMany(r => r).ToKMatrix(dim);
        list.Add((perm, cycles, mat));
    }

    return list.ToArray();
}

ConcreteGroup<KMatrix<ZnInt>> MetaCyclicGLnp_DiagByPerm_Slow(int m, int n, int r, int dim)
{
    var nks = Partitions32[dim].Select(l => l.Aggregate((a0, a1) => a0 * a1))
        .SelectMany(e => Dividors(e).Append(e).Where(j => j != 1))
        .Append(n)
        .Distinct()
        .ToArray();
    foreach (var p in nks.Select(nk => Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % nk == 0)).Order())
    {
        var Fp = FG.UnInt(p);
        var o = Fp.Neutral();
        var GL = FG.GLnK($"F{p}", dim, o);
        // Console.WriteLine($"Solve M({m}x:{n}){r} in {GL}");
        
        var m1s = Fp.Where(e => n % Fp.ElementsOrders[e] == 0)
            .SelectMany(e => PartitionsToMatrix(GL, e))
            .Where(m1 => IsOrder(m1.mat, n))
            .ToArray();

        foreach (var (perm, cycles, m1) in m1s)
        {
            var m1i = m1.Inv();
            var seq = cycles.Select(c => c.Length).Select(l =>
            {
                var r0 = IntExt.PowMod(r, l, m);
                var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0 && e.Pow(r0).Equals(e))
                    .OrderByDescending(e => Fp.ElementsOrders[e]).ToArray();

                return ordm.Select(a => l.Range(1).Select(k => a.Pow(IntExt.PowMod(r, k, m))).Reverse().ToArray());
            }).MultiLoop().Select(l => l.ToArray());
            foreach (var l in seq)
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

void Test()
{
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(4, 4, 3, 3));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(7, 3, 2, 3));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(5, 4, 2, 2));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(5, 4, 2, 3));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(5, 4, 2, 4));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(11, 10, 2, 5));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(11, 10, 2, 6));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(11, 10, 2, 7));
    DisplayGroup.HeadOrdersGenerators(MetaCyclicGLnp_DiagByPerm_Slow(11, 10, 2, 10));
}

void AllGensOfMtCycSdpUpToOrder(int maxOrd, int maxDim)
{
    GlobalStopWatch.Restart();
    Console.WriteLine("Start filtering MetaCyclic Groups");
    GlobalStopWatch.AddLap();
    var missing = new List<(int, int, int)>();
    var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => MetaCyclicSdp(ord))
        .Select(e => (e, FG.MetaCyclicSdp(e.m, e.n, e.r).AllSubgroups())).ToDictionary(e => e.Item2, e => e.e);
    var isoMtCycSdp = allMtCycSdp.Keys.FilterIsomorphic().ToDictionary(e => e, e => allMtCycSdp[e]);
    GlobalStopWatch.Show("End filtering");

    GlobalStopWatch.AddLap();
    foreach (var (m0, e) in isoMtCycSdp)
    {
        var id = FG.FindIdGroup(m0.Parent, m0.Infos);
        var found = false;
        for (int dim = 2; dim <= maxDim; dim++)
        {
            var mtGL = MetaCyclicGLnp_DiagByPerm_Slow(e.m, e.n, e.r, dim);
            if (mtGL.Count() != 1)
            {
                found = true;
                DisplayGroup.HeadOrdersGenerators(mtGL);
                if (!mtGL.IsIsomorphicTo(m0.Parent))
                    throw new();

                break;
            }
        }

        if (!found)
            missing.Add(e);
        else
        {
            if (id.Length != 0)
            {
                Console.WriteLine(id[0].FullName);
                Console.WriteLine();
            }
        }
    }

    var total = isoMtCycSdp.Count;
    missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    GlobalStopWatch.Show("End Gens");
    GlobalStopWatch.Show("END");
    Console.Beep();
}


{
    // AllGensOfMtCycSdpUpToOrder(maxOrd:64, maxDim:6);
    AllGensOfMtCycSdpUpToOrder(maxOrd:128, maxDim:10);
}