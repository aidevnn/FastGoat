using System.CodeDom;
using System.Collections;
using System.Diagnostics;
using System.Numerics;
using System.Reflection.Emit;
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
using FastGoat.UserGroup.DatabaseSmallGroups;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.GModuleN;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Words.Tools;
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

static MatFq[] DPGLnpGenerators(int n, int q)
{
    var og = FG.DPGLnpOrder(n, q);
    if (og > FG.MatrixGroupMaxOrder || PrimesDecomposition(q).Distinct().Count() != 1)
        throw new();

    var gl = new GLnq(n, q);
    Console.WriteLine(gl);
    var x = gl.Fq.X;
    var fq = Group.MulGroup("Fq", x);
    var e0 = fq.First(e => fq.ElementsOrders[e] == q - 1);
    var diag = gl.Neutral().Table.Select((c, i) => i == 0 ? e0 : c).ToArray();
    return FG.SnGensMat(n).Select(e => gl.Create(e.Table.Select(k => k * x.One).ToArray()))
        .Prepend(gl.Create(diag)).ToArray();
}

static MatFq[] DPSLnpGenerators(int n, int q)
{
    var og = FG.DPSLnpOrder(n, q);
    if (og > FG.MatrixGroupMaxOrder)
        throw new();

    var gl = new GLnq(n, q);
    var x = gl.Fq.X;
    var fq = Group.MulGroup("Fq", x);
    var e0 = fq.First(e => fq.ElementsOrders[e] == q - 1);
    var e1 = e0.Inv();
    var diag = gl.Neutral().Table.Select((c, i) => i == 0 ? e0 : (i == n + 1 ? e1 : c)).ToArray();
    return FG.SnGensMat(n).Select(e => gl.Create(e.Table.Select(k => k * x.One).ToArray())).Select(m =>
    {
        var deti = m.Det.Inv();
        var m0 = gl.Create(gl.Neutral().Table.Select((c, i) => i == 0 ? deti : c).ToArray());
        return gl.Op(m0, m);
    }).Prepend(gl.Create(diag)).ToArray();
}

List<int> FiniteFieldsOrder(int max)
{
    return Primes10000.Where(p => p <= max)
        .Select(p => (p, r: (int)(Math.Log10(max) / Math.Log10(p))))
        .SelectMany(e => e.r.Range(1).Select(k => e.p.Pow(k)))
        .Order()
        .ToList();
}

ConcreteGroup<MatFq> MetaCyclicGLnq_DiagByPerm(int m, int n, int r, int dim)
{
    var distinctTypes = IntExt.Partitions32[dim].Select(l => l.Order().ToArray()).OrderBy(l => l.Length).ToArray();
    var nks = distinctTypes.Select(l => l.Aggregate((a0, a1) => a0 * a1))
        .SelectMany(e => IntExt.Dividors(e).Append(e).Where(j => j != 1)).Append(n).ToHashSet();
    foreach (var q in nks.SelectMany(nk => FiniteFieldsOrder(1000).Where(q => (q - 1) % m == 0 && (q - 1) % nk == 0))
                 .Distinct().Order().Take(10))
    {
        var glq = new GLnq(dim, q);
        var Fq = Group.MulGroup($"F{q}", glq.Fq.X);
        var glqId = glq.Neutral().Table;
        var ordn = Fq.Where(e => n % Fq.ElementsOrders[e] == 0)
            .OrderBy(e => Fq.ElementsOrders[e])
            .Select(e => glqId.Select((c, i) => i == 0 ? e : c).ToArray())
            .Select(mat0 => glq.Create(mat0))
            .ToArray();
        var sn = new Sn(dim);
        var m1s = IntExt.Partitions32[dim].OrderBy(l => l.Count)
            .Select(t => IntExt.PermAndCyclesFromType(t.Order().ToArray()))
            .Select(e =>
            {
                var e0 = glq.Neutral().Table.Chunk(dim).ToArray();
                var perm = sn.CreateElement(e.perm.Select(i => i + 1).ToArray());
                var e1 = perm.Apply(e0);
                var mat0 = glq.Create(e1.SelectMany(v => v).ToArray());
                return ordn.Select(mat => glq.Op(mat0, mat))
                    .Where(mat => mat.IsOrder(n))
                    .Select(mat => (e.perm, e.cycles, mat));
            })
            .SelectMany(e => e);

        foreach (var (perm, cycles, m1) in m1s)
        {
            var m1i = glq.Invert(m1);
            var seq = cycles.Select(c => c.Length).Select(l =>
            {
                var r0 = IntExt.PowMod(r, l, m);
                var ordm = Fq.Where(e => m % Fq.ElementsOrders[e] == 0 && e.Pow(r0).Equals(e))
                    .OrderByDescending(e => Fq.ElementsOrders[e]);
                return ordm.Select(a => l.Range(1).Select(k => a.Pow(IntExt.PowMod(r, k, m))).Reverse().ToArray());
            }).MultiLoop().Select(l => l.ToArray());
            foreach (var l in seq)
            {
                var arr = (dim * dim).Range().Select(e => glq.Fq.X.Zero).ToArray();
                foreach (var (sols, idxs) in l.Zip(cycles))
                foreach (var (idx, sol) in idxs.Zip(sols))
                    arr[perm[idx] * (dim + 1)] = sol;

                var m0 = glq.Create(arr);
                if (m0.IsOrder(m) && glq.Op(m1i, glq.Op(m0, m1)).Equals(glq.Times(m0, r)))
                {
                    var name = IntExt.Gcd(m, n * (r - 1)) == 1 ? $"F({m}x:{n}){r}" : $"M({m}x:{n}){r}";
                    var mtGL = Group.Generate(name, glq, m0, m1);
                    if (mtGL.Count() == m * n)
                        return mtGL;
                }
            }
        }
    }

    return Group.Generate(new GLnq(1, 2));
}

ConcreteGroup<MatFq> MetaCyclicSdpMat(int m, int n, int r, int maxDim = 12)
{
    foreach (var dim in maxDim.Range(1).Where(d => d != 1))
    {
        var mtGL = MetaCyclicGLnq_DiagByPerm(m, n, r, dim);
        if (mtGL.Count() != 1)
            return mtGL;
    }

    return Group.Generate(new GLnq(1, 2));
}

void MtSdpGLnq()
{
    var ordMax = 64;
    DisplayGroup.HeadGenerators(FG.MetaCyclicSdpMat(8, 4, 7));
    foreach (var ord in ordMax.Range(2))
    {
        var seq = IntExt.Dividors(ord).Where(d => d > 1)
            .SelectMany(m => FG.MetaCyclicSdpGetR(m, ord / m).Select(r => (m, n: ord / m, r)))
            .Select(e => (e, MetaCyclicSdpMat(e.m, e.n, e.r, 5)))
            .Where(e => e.Item2.Count() != 1);
        foreach (var (e, mat1) in seq)
        {
            var mat2 = FG.MetaCyclicSdpMat(e.m, e.n, e.r);
            DisplayGroup.HeadGenerators(mat1);
            DisplayGroup.HeadGenerators(mat2);
            var n1q = (mat1.Neutral().GLnq.N, mat1.Neutral().GLnq.Fq.Q);
            var n2p = (mat2.Neutral().GL.N, mat2.Neutral().GL.P);
            if (n1q != n2p)
                Console.WriteLine($"#### Different n1q:{n1q} != n2p:{n2p}");
        }
    }
}

{
    // MtSdpGLnq();

    var seq = FiniteFieldsOrder(30);
    foreach (var dim in 4.Range(2))
    {
        foreach (var q in seq.Where(q => FG.DPGLnpOrder(dim, q) < 5000))
        {
            var gens = DPGLnpGenerators(dim, q);
            DisplayGroup.HeadGenerators(Group.Generate($"DPGL({dim},F{q})", gens[0].GLnq, gens));
        }
    }
    
    foreach (var dim in 4.Range(2))
    {
        foreach (var q in seq.Where(q => FG.DPSLnpOrder(dim, q) < 5000))
        {
            var gens = DPSLnpGenerators(dim, q);
            DisplayGroup.HeadGenerators(Group.Generate($"DPSL({dim},F{q})", gens[0].GLnq, gens));
        }
    }
}