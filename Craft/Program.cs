using System.Numerics;
using System.Text;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

ConcreteGroup<DPelt> MetaCyclicGLnp_DiagByPerm(int m, int n, int r, int dim)
{
    // Console.WriteLine($"Solve M({m}x:{n}){r} in GL({dim},{p})");
    var distinctTypes = IntExt.Partitions32[dim].Select(l => l.Order().ToArray()).OrderBy(l => l.Length).ToArray();
    var nks = distinctTypes.Select(l => l.Aggregate((a0, a1) => a0 * a1))
        .SelectMany(e => IntExt.Dividors(e).Append(e).Where(j => j != 1)).Append(n).ToHashSet();
    foreach (var p in nks.Select(nk => IntExt.Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % nk == 0))
                 .Distinct().Order())
    {
        var Up = FG.UnInt(p);
        var dpgl = new DPGL(dim, p);
        
        var ordn = Up.Where(e => n % Up.ElementsOrders[e] == 0)
            .OrderBy(e => Up.ElementsOrders[e])
            .Select(e => dpgl[dim.SeqLazy().Select(i => i == 0 ? e.K : 1).ToArray()])
            .ToArray();

        var m1s = IntExt.Partitions32[dim].OrderBy(l => l.Count)
            .Select(t => IntExt.PermAndCyclesFromType(t.Order().ToArray()))
            .Select<(int[] perm, int[][] cycles), IEnumerable<(int[] perm, int[][] cycles, DPelt mat)>>(e =>
            {
                var perm = dpgl.Sn.CreateElementTable(e.perm);
                return ordn.Select(dpelt => dpgl.Op(dpgl[perm], dpgl[dpelt.Diag]))
                    .Where(mat => mat.IsOrder(n))
                    .Select(mat => (e.perm, e.cycles, mat));
            })
            .SelectMany(e => e);

        foreach (var (perm, cycles, m1) in m1s)
        {
            var m1i = dpgl.Invert(m1);
            var seq = cycles.Select(c => c.Length).Select(l =>
            {
                var r0 = IntExt.PowMod(r, l, m);
                var ordm = Up.Where(e => m % Up.ElementsOrders[e] == 0 && e.Pow(r0).Equals(e))
                    .OrderByDescending(e => Up.ElementsOrders[e]);
                return ordm.Select(a => l.Range(1).Select(k => a.Pow(IntExt.PowMod(r, k, m))).Reverse().ToArray());
            }).MultiLoop().Select(l => l.ToArray());
            
            foreach (var l in seq)
            {
                var arr = new int[dim];
                foreach (var (sols, idxs) in l.Zip(cycles))
                foreach (var (idx, sol) in idxs.Zip(sols))
                    arr[perm[idx]] = sol.K;

                var m0 = dpgl[arr];
                if (m0.IsOrder(m) && dpgl.Op(m1i, dpgl.Op(m0, m1)).Equals(dpgl.Times(m0, r)))
                {
                    var name = IntExt.Gcd(m, n * (r - 1)) == 1 ? $"F({m}x:{n}){r}" : $"M({m}x:{n}){r}";
                    var mtGL = Group.Generate(name, dpgl, m0, m1);
                    if (mtGL.Count() == m * n)
                        return mtGL;
                }
            }
        }
    }

    return Group.Generate(new DPGL(1, 2));
}

ConcreteGroup<DPelt> MetaCyclicSdpMat(int m, int n, int r, int maxDim = 12)
{
    foreach (var dim in maxDim.Range(1).Where(d => d != 1 && (IntExt.Gcd(m, d) != 1 || IntExt.Gcd(m - 1, d) != 1)))
    {
        var mtGL = MetaCyclicGLnp_DiagByPerm(m, n, r, dim);
        if (mtGL.Count() != 1)
            return mtGL;
    }

    throw new GroupException(GroupExceptionType.GroupDef);
}

void TestMetaCyclic(int ord)
{
    foreach (var (m, n, r) in IntExt.Dividors(ord).Where(d => d > 1)
                 .SelectMany(m => FG.MetaCyclicSdpGetR(m, ord / m).Select(r => (m, n: ord / m, r))))
    {
        var mat1 = MetaCyclicSdpMat(m, n, r);
        DisplayGroup.HeadGenerators(mat1);
        var mat2 = FG.MetaCyclicSdpMat(m, n, r);
        DisplayGroup.HeadGenerators(mat2);
        if (!mat2.IsIsomorphicTo(mat1))
        {
            DisplayGroup.HeadNames(mat1);
            DisplayGroup.HeadNames(mat2);
            throw new();
        }
    }
}

DPelt[] DPGLnpGenerators(int n, int p)
{
    var og = FG.DPGLnpOrder(n, p);
    if (og > FG.MatrixGroupMaxOrder)
        throw new();

    if (!IntExt.Primes10000.Contains(p))
        throw new($"p={p} isnt prime");
    
    var dpgl = new DPGL(n, p);
    var e0 = NumberTheory.PrimitiveRootMod(dpgl.P);
    var id = Enumerable.Repeat(1, n).ToArray();
    var diag = dpgl[n.SeqLazy().Select(i => i == 0 ? e0 : 1).ToArray()];
    return dpgl.Sn.GetGenerators().Select(e => dpgl[id, e]).Prepend(diag).ToArray();
}

DPelt[] DPSLnpGenerators(int n, int p)
{
    var og = FG.DPSLnpOrder(n, p);
    if (og > FG.MatrixGroupMaxOrder)
        throw new($"og = {og}");

    if (!IntExt.Primes10000.Contains(p))
        throw new($"p={p} isnt prime");

    var dpgl = new DPGL(n, p);
    var gl = dpgl.GL;
    var e0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(gl.P, p - 1);
    var e1 = IntExt.InvModPbez(e0, p);
    var diag = dpgl[n.SeqLazy().Select(i => i == 0 ? e0 : i == 1 ? e1 : 1).ToArray()];
    return dpgl.Sn.GetGenerators().Select(perm =>
    {
        var det = MatrixExt.DeterminantByPivot(n, p, MatrixExt.Permutation(perm.Table));
        var dg = n.SeqLazy().Select(i => i == 0 ? IntExt.InvModPbez(det, p) : 1).ToArray();
        return dpgl[dg, perm];
    }).Prepend(diag).ToArray();
}

ConcreteGroup<DPelt> DPSLnp(int n, int p)
{
    var gens = DPSLnpGenerators(n, p);
    var dpgl = gens[0].Dpgl;
    return Group.Generate($"DPSL({n},{p})", dpgl, gens);
}

ConcreteGroup<DPelt> DPGLnp(int n, int p)
{
    var gens = DPGLnpGenerators(n, p);
    var dpgl = gens[0].Dpgl;
    return Group.Generate($"DPGL({n},{p})", dpgl, gens);
}

void Dpgl()
{
    var (n, p) = (6, 3);
    var dpgl1 = FG.DPGLnp(n, p);
    DisplayGroup.HeadGenerators(dpgl1);
    
    var dpgl2 = DPGLnp(n, p);
    DisplayGroup.HeadGenerators(dpgl2);

    DisplayGroup.AreIsomorphics(dpgl1, dpgl2);

    GlobalStopWatch.Bench(5, "MatMul", () => FG.DPGLnp(n, p));
    GlobalStopWatch.Bench(5, "MonomialMat", () => DPGLnp(n, p));
    GlobalStopWatch.Bench(5, "MatMul", () => FG.DPGLnp(n, p));
    GlobalStopWatch.Bench(5, "MonomialMat", () => DPGLnp(n, p));
}

void Dpsl()
{
    var (n, p) = (6, 3);
    var dpsl1 = FG.DPSLnp(n, p);
    DisplayGroup.HeadGenerators(dpsl1);
    
    var dpsl2 = DPSLnp(n, p);
    DisplayGroup.HeadGenerators(dpsl2);

    // DisplayGroup.AreIsomorphics(dpsl1, dpsl2);

    GlobalStopWatch.Bench(5, "MatMul", () => FG.DPSLnp(n, p));
    GlobalStopWatch.Bench(5, "MonomialMat", () => DPSLnp(n, p));
    GlobalStopWatch.Bench(5, "MatMul", () => FG.DPSLnp(n, p));
    GlobalStopWatch.Bench(5, "MonomialMat", () => DPSLnp(n, p));
}

void Bench(int m, int n, int r)
{
    var mat1 = MetaCyclicSdpMat(m, n, r);
    DisplayGroup.HeadGenerators(mat1);
    foreach (var mat in mat1.GetGenerators())
        mat.DisplayArrays();
    
    var mat2 = FG.MetaCyclicSdpMat(m, n, r);
    DisplayGroup.HeadGenerators(mat2);
    DisplayGroup.AreIsomorphics(mat1, mat2);
    
    GlobalStopWatch.Bench(5, "MatMul", () => FG.MetaCyclicSdpMat(m, n, r));
    GlobalStopWatch.Bench(5, "MonomialMat", () => MetaCyclicSdpMat(m, n, r));
    GlobalStopWatch.Bench(5, "MatMul", () => FG.MetaCyclicSdpMat(m, n, r));
    GlobalStopWatch.Bench(5, "MonomialMat", () => MetaCyclicSdpMat(m, n, r));
}

void TestOrd157()
{
    var dpgl = new DPGL(12, 53);
    var m1a = dpgl[[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1)];
    var m1b = m1a.ToGL();
    var m2a = dpgl[[44, 16, 49, 46, 24, 36, 47, 10, 13, 15, 42, 28], (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)];
    var m2b = m2a.ToGL();
    var mat1 = Group.Generate("F(13x:12)2", dpgl, m1a, m2a);
    DisplayGroup.HeadGenerators(mat1);
    var mat2 = Group.Generate("F(13x:12)2", dpgl.GL, m1b, m2b);
    DisplayGroup.HeadGenerators(mat2);
    
    GlobalStopWatch.Bench(5, "MonomialMat", () => Group.Generate("F(13x:12)2", dpgl, m1a, m2a));
    GlobalStopWatch.Bench(5, "MatMul", () => Group.Generate("F(13x:12)2", dpgl.GL, m1b, m2b));
    GlobalStopWatch.Bench(5, "MonomialMat", () => Group.Generate("F(13x:12)2", dpgl, m1a, m2a));
    GlobalStopWatch.Bench(5, "MatMul", () => Group.Generate("F(13x:12)2", dpgl.GL, m1b, m2b));
    GlobalStopWatch.Bench(50, "MonomialMat", () => Group.Generate("F(13x:12)2", dpgl, m1a, m2a));
    GlobalStopWatch.Bench(50, "MatMul", () => Group.Generate("F(13x:12)2", dpgl.GL, m1b, m2b));
}

{
    // TestOrd157();
    // Bench(13, 12, 2);
    //
    // for (int k = 20; k < 128; ++k)
    //     TestMetaCyclic(k);
    //
    // Dpgl();
    Dpsl();

}