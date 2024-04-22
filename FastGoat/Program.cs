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
Ring.DisplayPolynomial = MonomDisplay.StarPowFct;

IEnumerable<Mat> GL2Candidate(GL gl, ZnInt e, int ord) => new[]
    {
        gl[0, e.K, 1, 0]
    }
    .Where(mat => IsOrder(mat, ord));

IEnumerable<Mat> GL3Candidate(GL gl, ZnInt e, int ord) => new[]
    {
        gl[e.K, 0, 0, 0, 0, 1, 0, 1, 0],
        gl[0, e.K, 0, 0, 0, 1, 1, 0, 0]
    }
    .Where(mat => IsOrder(mat, ord));

IEnumerable<Mat> GL4Candidate(GL gl, ZnInt e, int ord) => new[]
    {
        gl[0, e.K, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0],
        gl[e.K, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0]
    }
    .Where(mat => IsOrder(mat, ord));

IEnumerable<Mat> GL5Candidate(GL gl, ZnInt e, int ord) => new[]
    {
        gl[0, e.K, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
        gl[e.K, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
    }
    .Where(mat => IsOrder(mat, ord));

IEnumerable<Mat> GL6Candidate(GL gl, ZnInt e, int ord) => new[]
    {
        gl[0, e.K, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
        gl[0, e.K, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
        gl[0, e.K, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
    }
    .Where(mat => IsOrder(mat, ord));

IEnumerable<Mat> GLnCandidates(GL gl, ZnInt e, int ord)
{
    if (gl.N == 2)
        return GL2Candidate(gl, e, ord);
    if (gl.N == 3)
        return GL3Candidate(gl, e, ord);
    if (gl.N == 4)
        return GL4Candidate(gl, e, ord);
    if (gl.N == 5)
        return GL5Candidate(gl, e, ord);
    if (gl.N == 6)
        return GL6Candidate(gl, e, ord);

    throw new();
}

IEnumerable<Mat> Candidat(int dim, int n, int p)
{
    var Fp = FG.UnInt(p);
    var gl = new GL(dim, p);
    var ordn = Fp.Where(e => n % Fp.ElementsOrders[e] == 0).ToArray();
    return ordn.SelectMany(e => GLnCandidates(gl, e, n));
}

Dictionary<int, ZnInt> SolveSys(int m, Polynomial<ZnInt, Xi>[] Sys, Dictionary<int, ZnInt> SolsF, int depth = 0)
{
    if (depth > 6)
    {
        Sys.Println();
        SolsF.Println();
        throw new();
    }

    if (Sys.Length == 0)
        return SolsF;

    var z = Sys[0].Zero;
    var p = z.P;
    var Fp = FG.UnInt(p);

    var oneUnknown = Sys.Where(P => P.NbIndeterminates == 1).ToArray();
    var ind = z.Indeterminates.ToList();
    var sols = new Dictionary<Xi, Polynomial<ZnInt, Xi>[]>();
    var solsF = new Dictionary<int, ZnInt>(SolsF);
    var sys0 = Sys.ToList();
    foreach (var P in oneUnknown)
    {
        var xi = P.ExtractIndeterminate;
        sols[xi] = Fp.Where(c => m % Fp.ElementsOrders[c] == 0)
            .Select(c => c * z.One)
            .Where(c => P.Substitute(c, xi).IsZero())
            .OrderByDescending(c => Fp.ElementsOrders[c.ConstTerm])
            .ThenBy(c => c)
            .ToArray();
        solsF[ind.FindIndex(xj => xj.Equals(xi))] = sols[xi][0].ConstTerm;
        sys0.Remove(P);
    }

    var sys1 = sys0.Select(P => P.Substitute(sols.Select(e => (e.Key, e.Value[0])).ToList()))
        .Where(P => !P.IsZero())
        .ToList();

    // sols.Select(e => $"{e.Key} in [{e.Value.Glue(", ")}]").Println("One Unknown");
    // sys1.Println("Remaining system");
    if (sys1.Count == 1 && sys1[0].Degree == 0)
        return z.Indeterminates.Length.Range().ToDictionary(i => i, i => Fp.Neutral());

    var xis = sys1.SelectMany(P0 => P0.ExtractAllIndeterminates).ToArray();
    var xis0 = xis.OrderByDescending(xi => sys1.Count(P0 => P0.ExtractAllIndeterminates.Contains(xi))).ToArray();
    if (xis0.Length != 0)
    {
        var xi = xis0[0];
        var c0 = Fp.First(e => Fp.ElementsOrders[e] == m);
        sols[xi] = [c0 * z.One];
        solsF[ind.FindIndex(xj => xj.Equals(xi))] = sols[xi][0].ConstTerm;
        sys1 = sys0.Select(P => P.Substitute(sols.Select(e => (e.Key, e.Value[0])).ToList()))
            .Where(P => !P.IsZero())
            .ToList();

        // sys1.Println($"Remaining system with {xi} == {c0}");
    }

    return SolveSys(m, sys1.ToArray(), solsF, depth + 1);
}

ConcreteGroup<Mat> SymbolicSolveMetaCyclic(int m, int n, int r, int dim)
{
    var nc = n % dim != 0 ? n : n / dim;
    var p = Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % nc == 0);
    return SymbolicSolveMetaCyclic2(m, n, r, dim, p);
}

ConcreteGroup<Mat> SymbolicSolveMetaCyclic2(int m, int n, int r, int dim, int p)
{
    // Console.WriteLine($"Solve M({m}x:{n}){r} Symbolic in GL({dim},{p})");
    var Fp = FG.UnInt(p);

    var gen = Fp.GetGenerators().First();
    var xs = Ring.Polynomial(gen, MonomOrder.Lex, (dim, "x")).Deconstruct();
    var bs = new PolynomialBasis<ZnInt, Xi>(xs.SelectMany(P => new[] { P.Pow(m) - 1, P.Pow(p) - 1 }).ToArray());
    var Xs = xs.Select(P => new EPolynomial<ZnInt>(P, P.One, bs)).Deconstruct();
    var (z, o) = (Xs[0].Zero, Xs[0].One);

    var M0 = dim.Range().Grid2D().Select(e => e.t1 == e.t2 ? Xs[e.t1] : z).ToKMatrix(dim);
    var M0_r = M0.Pow(r);

    foreach (var m1 in Candidat(dim, n, p))
    {
        var M1 = m1.Table.Select(c => c * o).ToKMatrix(dim);
        if (!M1.Pow(n).Equals(M1.One))
            throw new("Problem1");

        var Sys = (M1.Inv() * M0 * M1 - M0_r).Where(P => !P.IsZero()).Select(P => P.Num).ToArray();
        var Sols = Ring.ReducedGrobnerBasis(Sys);
        if (Sols.Where(P => !P.IsZero()).All(P => P.Degree > 1))
        {
            var Sols0 = Sols.Select(P => new EPolynomial<ZnInt>(P, P.One, bs)).Where(P => !P.IsZero()).Select(P => P.Num).ToArray();
            var solsF = SolveSys(m, Sols0, new());

            // var ind = Sols0[0].Indeterminates;
            // Console.WriteLine("M0");
            // Console.WriteLine(M0);
            // Console.WriteLine($"M0^{m} = I{dim} : {M0.Pow(m).Equals(M0.One)}");
            // Console.WriteLine();
            // Console.WriteLine("M1");
            // Console.WriteLine(M1);
            // Console.WriteLine($"M1^{n} = I{dim} : {M1.Pow(n).Equals(M1.One)}");
            // Console.WriteLine();
            // Sys.Println($"Sys in GL({dim}, {p}), M0^{m}=Id, M1^{n}=Id, and M1^-1 * M0 * M1 = M1^{r}");
            // Sols0.Println("Sols");
            // Console.WriteLine();
            // solsF.Select(e => $"{ind[e.Key]} = {e.Value}").Println("Final Solutions");

            var arr = new int[dim * dim];
            for (int i = 0; i < dim; i++)
                arr[i * dim + i] = solsF[i].K;

            var m0 = m1.GL.Create(arr);
            if (IsOrder(m0, m))
            {
                var mtGL = Group.Generate($"M({m}x:{n}){r}", m0.GL, m0, m1);
                if (mtGL.Count() == m * n)
                    return mtGL;
            }
            else
            {
                // Console.WriteLine("M0 order problem");
            }
        }
        else
        {
            // Console.WriteLine("No M0 found");
        }
    }

    return Group.Generate(new GL(1, 2));
}

(int m, int n, int r)[] MetaCyclicSdp(int order)
{
    return IntExt.Dividors(order).Where(d => d > 1)
        .SelectMany(m => FG.MetaCyclicSdpGetR(m, order / m).Select(r => (m, n: order / m, r)))
        .ToArray();
}

(int m, int n, int r)[] FrobeniusSdp(int order)
{
    return IntExt.Dividors(order).Where(d => d > 1)
        .SelectMany(m => FG.FrobeniusGetR(m, order / m).Select(r => (m, n: order / m, r)))
        .ToArray();
}

bool IsOrder(Mat m, int o)
{
    var gl = m.GL;
    var e = gl.Neutral();
    var mk = gl.Neutral();
    if (PowMod(m.Det, o, gl.P) != 1)
        return false;

    for (int k = 0; k < o; k++)
    {
        if (k != 0 && mk.Equals(e))
            return false;

        mk = gl.Op(mk, m);
    }

    return mk.Equals(e);
}

ConcreteGroup<Mat> MetaCyclicGL2p_Meth1(int m, int n, int r)
{
    var mtSdp = FG.MetaCyclicSdp(m, n, r);
    var p = Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % n == 0);

    var gl = new GL(2, p);
    var Fp = FG.UnInt(p);
    var ordms = Fp.Where(e => Fp.ElementsOrders[e] == m).Select(e => e.K).ToArray();
    var ordns = Fp.Where(e => Fp.ElementsOrders[e] == n).Select(e => e.K).ToArray();

    var a0s = ordms.Grid2D().Where(e => PowMod(e.t1, r, p) == e.t2 && PowMod(e.t2, r, p) == e.t1).ToArray();
    var b0s = ordns.Where(e => IsOrder(gl[0, e, 1, 0], n)).ToArray();

    foreach (var (a0, a1) in a0s)
    {
        foreach (var b0 in b0s)
        {
            {
                var m0 = gl[a0, 0, 0, a1];
                var m1 = gl[0, b0, 1, 0];
                var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
                if (mtGl.Count() == m * n)
                    return mtGl;

                // if (mtGl.IsIsomorphicTo(mtSdp))
                //     return mtGl;
            }
        }
    }

    return Group.Generate(new GL(1, 2));
}

ConcreteGroup<Mat> MetaCyclicGL2p_Meth2(int m, int n, int r)
{
    var mtSdp = FG.MetaCyclicSdp(m, n, r);
    foreach (var p in Primes10000.Where(p => m % p == 0 && (p - 1) % n == 0))
    {
        var Fp = FG.UnInt(p);
        var gl = new GL(2, p);

        var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0).ToArray();
        var ordn = Fp.Where(e => Fp.ElementsOrders[e] == n).ToArray();

        var m0s = ordm.Grid3D(ordm, Fp.Prepend(Fp.Neutral().Zero).OrderBy(e => e.K))
            .Select(e => gl[e.t1.K, e.t3.K, 0, e.t2.K])
            .Where(mat => IsOrder(mat, m))
            .ToArray();

        var m1s = ordn.Grid2D(ordn.Append(Fp.Neutral())).Select(e => gl[e.t1.K, 0, 0, e.t2.K])
            .Where(mat => IsOrder(mat, n))
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

IEnumerable<Mat> GL3MatrixPerm(GL gl, (ZnInt t1, ZnInt t2, ZnInt t3) e, int ord) => new[]
    {
        gl[e.t1.K, 0, 0, 0, 0, e.t2.K, 0, e.t3.K, 0],
        gl[0, e.t1.K, 0, 0, 0, e.t2.K, e.t3.K, 0, 0]
    }
    .Where(mat => IsOrder(mat, ord));

ConcreteGroup<Mat> MetaCyclicGL3p_Meth(int m, int n, int r)
{
    var mtSdp = FG.MetaCyclicSdp(m, n, r);
    var p = Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % n == 0);
    var gl = new GL(3, p);
    var Fp = FG.UnInt(p);
    var ordm = Fp.Where(e => m % Fp.ElementsOrders[e] == 0).OrderBy(e => e.K).ToArray();
    var ordn = Fp.Where(e => n % Fp.ElementsOrders[e] == 0).OrderBy(e => e.K).ToArray();

    var m0s = ordm.Grid3D().Select(e => gl[e.t1.K, 0, 0, 0, e.t2.K, 0, 0, 0, e.t3.K])
        .Where(mat => IsOrder(mat, m))
        .ToArray();
    var m1s = ordn.Grid3D().SelectMany(e => GL3MatrixPerm(gl, e, n)).ToArray();

    foreach (var (m0, m1) in m0s.Grid2D(m1s).Where(e => gl.Op(gl.Invert(e.t2), gl.Op(e.t1, e.t2)).Equals(gl.Times(e.t1, r))))
    {
        var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
        if (mtGl.Count() == m * n)
            return mtGl;
    }

    return Group.Generate(new GL(1, 2));
}

ConcreteGroup<Mat> MetaCyclicGL4p_Meth(int m, int n, int r)
{
    if (n % 4 != 0)
        return Group.Generate(new GL(1, 2));

    var p = Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % (n / 4) == 0);

    var Fp = FG.UnInt(p);
    var Fp4 = Fp.Where(e => e.Pow(m).K == 1)
        .MultiLoop(4)
        .Select(l => l.ToArray()).Select(l => (t1: l[0], t2: l[1], t3: l[2], t4: l[3]))
        .OrderBy(l => l)
        .ToArray();

    var gl = new GL(4, p);
    var m1s = Fp.Where(e => Fp.ElementsOrders[e] == n / 4)
        .SelectMany(e => new[]
        {
            gl[0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, e.K, 0],
            gl[0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, e.K, 0, 0, 0]
        })
        .Where(mat => IsOrder(mat, n))
        .ToArray();

    var m0s = Fp4.Select(e => gl[e.t1.K, 0, 0, 0, 0, e.t2.K, 0, 0, 0, 0, e.t3.K, 0, 0, 0, 0, e.t4.K])
        .Where(mat => IsOrder(mat, m))
        .ToArray();

    foreach (var (m0, m1) in m0s.Grid2D(m1s).Where(e => gl.Op(gl.Invert(e.t2), gl.Op(e.t1, e.t2)).Equals(gl.Times(e.t1, r))))
    {
        var mtGl = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
        if (mtGl.Count() == m * n)
            return mtGl;
    }

    return Group.Generate(new GL(1, 2));
}

void Run(int maxOrd = 32, bool frob = false)
{
    GlobalStopWatch.Restart();
    Func<int, (int m, int n, int r)[]> f = frob ? FrobeniusSdp : MetaCyclicSdp;

    var missing = new List<(int, int, int)>();
    var allMtCycSdp = (maxOrd - 5).Range(6).SelectMany(ord => f(ord))
        .Select(e => (e, FG.MetaCyclicSdp(e.m, e.n, e.r).AllSubgroups())).ToDictionary(e => e.Item2, e => e.e);
    var isoMtCycSdp = allMtCycSdp.Keys.FilterIsomorphic().ToDictionary(e => e, e => allMtCycSdp[e]);

    foreach (var (m0, e) in isoMtCycSdp)
    {
        var id = FG.FindIdGroup(m0.Parent, m0.Infos);
        var mtGL1 = MetaCyclicGL2p_Meth1(e.m, e.n, e.r);
        if (mtGL1.Count() != 1)
        {
            DisplayGroup.HeadOrdersGenerators(mtGL1);
            if (id.Length != 0)
            {
                Console.WriteLine(id[0].FullName);
                Console.WriteLine();
            }
            continue;
        }

        var mtGL0 = MetaCyclicGL2p_Meth2(e.m, e.n, e.r);
        if (mtGL0.Count() != 1)
        {
            DisplayGroup.HeadOrdersGenerators(mtGL0);
            if (id.Length != 0)
            {
                Console.WriteLine(id[0].FullName);
                Console.WriteLine();
            }
            continue;
        }

        var found = false;
        for (int dim = 2; dim <= 6; dim++)
        {
            var mtGL = SymbolicSolveMetaCyclic(e.m, e.n, e.r, dim);
            if (mtGL.Count() != 1)
            {
                found = true;
                DisplayGroup.HeadOrdersGenerators(mtGL);
                if (id.Length != 0)
                {
                    Console.WriteLine(id[0].FullName);
                    Console.WriteLine();
                }
                break;
            }
        }

        if (!found)
            missing.Add(e);
    }

    var total = isoMtCycSdp.Count;
    missing.Println(e => $"M({e.Item1}x:{e.Item2}){e.Item3}", $"Missing:{missing.Count} Found:{total - missing.Count}/{total}");
    Console.WriteLine($"var missing = new [] {{ {missing.Glue(", ")} }};");
    GlobalStopWatch.Show("END");
    Console.Beep();
}

{
    Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracketNoFmt;
    Group.ActivedStorage(false);

    // Run(frob: true);
    //
    // Run(maxOrd:64);
    // Run(maxOrd:96);
    Run(128);
    // Run(128);
}
