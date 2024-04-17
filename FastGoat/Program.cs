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

(int m, int n, int r)[] MetaCyclicSdp(int order)
{
    var ms = IntExt.Dividors(order).Where(d => d > 1).ToArray();

    HashSet<(int m, int n, int r)> all = new();
    foreach (var m in ms)
    {
        var n = order / m;
        var rs = IntExt.SolveAll_k_pow_m_equal_one_mod_n(m, n);
        foreach (var r in rs)
            all.Add((m, n, r));
    }

    return all.ToArray();
}

(int, int, int, int) MatOp(int p, (int a00, int a01, int a10, int a11) e0, (int a00, int a01, int a10, int a11) e1)
{
    return ((e0.a00 * e1.a00 + e0.a01 * e1.a10) % p,
        (e0.a00 * e1.a01 + e0.a01 * e1.a11) % p,
        (e0.a10 * e1.a00 + e0.a11 * e1.a10) % p,
        (e0.a10 * e1.a01 + e0.a11 * e1.a11) % p);
}

int OrdMat(int p, (int a00, int a01, int a10, int a11) e0, int max = 32)
{
    var e1 = (1, 0, 0, 1);
    var tmp = e0;
    var ord = 1;
    while (tmp != e1)
    {
        ord++;
        if (ord > max)
            return -1;
        tmp = MatOp(p, tmp, e0);
    }

    return ord;
}

Dictionary<int, LinkedList<(int, int, int, int)>> FastGL2p(int p)
{
    var a0s = SolveAll_k_pow_m_equal_one_mod_n(p, p - 1);
    var rg = p.Range().Grid3D(a0s, a0s).ToArray();
    var gl = new Dictionary<int, LinkedList<(int, int, int, int)>>();
    var c4 = FG.PermGroup(4, (1, 2, 3, 4));
    foreach (var (e0, e1, e2, e3) in rg.SelectMany(e => c4.Select(perm => perm.Apply(new[] { e.t1, e.t2, e.t3, 0 })))
                 .Select(f => (f[0], f[1], f[2], f[3])))
    {
        // if (e0 * e1 * e2 * e3 != 0)
        //     continue;

        var e = (e0, e1, e2, e3);
        var det = AmodP(e0 * e3 - e1 * e2, p);
        // Console.WriteLine($"{e} det {det}");
        if (det == 0)
            continue;

        var ord = OrdMat(p, e);
        if (ord == -1)
            continue;

        if (!gl.ContainsKey(ord))
            gl[ord] = new();

        gl[ord].AddLast(e);
    }

    return gl;
}

Dictionary<int, LinkedList<(int, int, int, int)>> ShowGLorders(int p)
{
    var gl = FastGL2p(p);
    var ord = gl.Sum(e => e.Value.Count);
    var ord2 = FG.GLnqOrder(2, p);
    Console.WriteLine($"|GL(2,{p})| = {ord}");
    // if (ord != ord2)
    //     throw new($"actual:{ord} expected:{ord2}");

    Console.WriteLine(
        $"Elements Orders : {gl.ToDictionary(e => e.Key, e => e.Value.Count).AscendingByKey().GlueMap(fmt: "[{0}] = {1}")}");
    Console.WriteLine();
    return gl;
}

Mat MatCreate(GL gl, (int a00, int a01, int a10, int a11) e0) => gl.Create(e0.a00, e0.a01, e0.a10, e0.a11);

bool IsOrder(Mat mat, int ord)
{
    var gl = mat.GL;
    var e = gl.Neutral();
    var tmp = gl.Neutral();
    for (int k = 0; k < ord; k++)
    {
        if (k != 0 && tmp.Equals(e))
            return false;

        tmp = gl.Op(tmp, mat);
    }

    return tmp.Equals(e);
}

(int m, int n, int r)[] FrobeniusSdp(int order)
{
    var ms = IntExt.Dividors(order).Where(d => d > 1 && d % 2 == 1).ToArray();

    List<(int m, int n, int r)> all = new();
    foreach (var m in ms)
    {
        var n = order / m;
        foreach (var r in FG.FrobeniusGetR(m, n))
            all.Add((m, n, r));
    }

    return all.ToArray();
}

bool BruteForce1(ConcreteGroup<Mat> gl, int m, int n, int r)
{
    var p = gl.Neutral().GL.P;
    var a0 = Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
    var a = new ZnInt(p, a0);
    var Zp = Group.MulGroup($"F{p}", a);
    var m0s = gl.Where(e => gl.ElementsOrders[e] == m).ToArray();
    var m1s = gl.Where(e => gl.ElementsOrders[e] == n && e.IsSym).ToArray();
    Console.WriteLine($"###### Start search generators of M({m}x:{n}){r}");
    Console.WriteLine($"     in {gl} Possibilities:{m0s.Length * m1s.Length}");
    foreach (var m0 in m0s)
    {
        var m0r = gl.Times(m0, r);
        if (m0r.Det != m0.Det)
            continue;
        foreach (var m1 in m1s)
        {
            var m1i = gl.Invert(m1);
            if (!gl.Op(m1i, gl.Op(m0, m1)).Equals(m0r))
                continue;

            var gmn = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
            if (gmn.Count() != m * n)
                continue;

            // Console.WriteLine($"###### Start search generators of M({m}x:{n}){r}");
            // Console.WriteLine($"     in {gl} Possibilities:{m0s.Length * m1s.Length}");
            DisplayGroup.HeadGenerators(gmn);
            Console.WriteLine();
            return true;
        }
    }

    return false;
}

bool BruteForce2(int p, LinkedList<(int a00, int a01, int a10, int a11)> ord_m, LinkedList<(int a00, int a01, int a10, int a11)> ord_n,
    int m, int n, int r)
{
    var gl = new GL(2, p);
    var tot_n = ord_n.Count * ord_m.Count;
    var m0s = ord_m.Select(e => MatCreate(gl, e)).ToArray();
    var m1s = ord_n.Select(e => MatCreate(gl, e)).ToArray();
    Console.WriteLine($"###### Start search generators of M({m}x:{n}){r}");
    Console.WriteLine($"     in {gl} Possibilities:{m0s.Length * m1s.Length}");
    foreach (var m0 in m0s)
    {
        var m0r = gl.Times(m0, r);
        if (m0r.Det != m0.Det)
            continue;

        foreach (var m1 in m1s)
        {
            var m1i = gl.Invert(m1);
            if (!gl.Op(m1i, gl.Op(m0, m1)).Equals(m0r))
                continue;

            var gmn = Group.Generate($"M({m}x:{n}){r}", gl, m0, m1);
            if (gmn.Count() != m * n)
                continue;

            // Console.WriteLine($"###### Start search generators of M({m}x:{n}){r}");
            // Console.WriteLine($"     in {gl} Possibilities:{m0s.Length * m1s.Length}");
            DisplayGroup.HeadGenerators(gmn);
            Console.WriteLine();
            return true;
        }
    }

    return false;
}

void Run1()
{
    GlobalStopWatch.AddLap();
    var gls0 = new List<ConcreteGroup<Mat>>();
    foreach (var p in Primes10000.Take(8))
    {
        var (a, b) = FG.GLnpGenerators(2, p);
        var gl = Group.Generate($"GL(2,{p})", a.GL, a, b);
        gls0.Add(gl);
        DisplayGroup.HeadOrders(gl);
    }

    {
        var (a, b) = FG.GLnpGenerators(4, 2);
        var gl = Group.Generate(a.GL, a, b);
        gls0.Add(gl);
        DisplayGroup.HeadOrders(gl);
    }

    GlobalStopWatch.Show();

    GlobalStopWatch.AddLap();
    var (nbFound, nbTotal) = (0, 0);
    var missing = new List<(int m, int n, int r)>();

    for (int ord = 6; ord < 65; ord++)
    {
        foreach (var (m, n, r) in MetaCyclicSdp(ord))
        {
            var found = false;
            if ((n == 2 && r == m - 1) || (n == 4 && m % 2 == 1 && r == m - 1))
                continue;
            else
            {
                ++nbTotal;
                foreach (var gl in gls0)
                {
                    if (BruteForce1(gl, m, n, r))
                    {
                        found = true;
                        ++nbFound;
                        break;
                    }
                }
            }

            if (!found)
            {
                missing.Add((m, n, r));
                // Console.WriteLine($"########################### generators of M({m}x:{n}){r} not found ###########################\n");
            }
        }
    }

    GlobalStopWatch.Show($"Found:{nbFound} Total:{nbTotal}");
    GlobalStopWatch.Show("End");
    missing.Println(e => $"M({e.m}x:{e.n}){e.r}", $"Count:{missing.Count}");
    Console.Beep();
}

void Run2()
{
    GlobalStopWatch.Restart();
    GlobalStopWatch.AddLap();
    var nbPrimes = 20;
    var dicoOrdsGl = new Dictionary<int, Dictionary<int, LinkedList<(int, int, int, int)>>>(nbPrimes);
    foreach (var p in Primes10000.Take(nbPrimes))
        dicoOrdsGl[p] = ShowGLorders(p);

    GlobalStopWatch.Show();

    GlobalStopWatch.AddLap();
    var (nbFound, nbTotal) = (0, 0);
    var missing = new List<(int m, int n, int r)>();

    for (int ord = 6; ord < 65; ord++)
    {
        foreach (var (m, n, r) in MetaCyclicSdp(ord))
        {
            var found = false;
            if ((n == 2 && r == m - 1) || (n == 4 && m % 2 == 1 && r == m - 1))
                continue;
            else
            {
                ++nbTotal;
                foreach (var p in Primes10000.Take(nbPrimes))
                {
                    if (!dicoOrdsGl[p].ContainsKey(m) || !dicoOrdsGl[p].ContainsKey(n))
                        continue;

                    if (BruteForce2(p, dicoOrdsGl[p][m], dicoOrdsGl[p][n], m, n, r))
                    {
                        found = true;
                        ++nbFound;
                        break;
                    }
                }
            }

            if (!found)
            {
                missing.Add((m, n, r));
                // Console.WriteLine($"########################### generators of M({m}x:{n}){r} not found ###########################\n");
            }
        }
    }

    GlobalStopWatch.Show($"Found:{nbFound} Total:{nbTotal}");
    GlobalStopWatch.Show("End");
    missing.Println(e => $"M({e.m}x:{e.n}){e.r}", $"Count:{missing.Count}");
    Console.Beep();
}

{
    // Run2();
    // Run2();
    // Run1();
}

void FrobPerm(int m, int n, int r)
{
    var mx = int.Max(m, n);

    foreach (var k in new[] { 1, int.Min(m, n) })
    {
        var sn = new Sn(k * mx);
        var km = sn.N / m;

        var cm = sn.ComposesCycles(km.Range().Select(i => new Tuple2Array(m.Range(1 + i * m))).ToArray());
        var cmr = sn.Times(cm, r);

        var N = new Dictionary<int, int>();
        var infos = SearchStep(n, cm, cmr, N);
        if (infos.state == 1)
        {
            var cn = sn.Invert(sn.CreateElement(sn.N.Range().Select(i => infos.cycle[i] + 1).ToArray()));
            var G = Group.Generate("G", sn, cm, cn);
            var sg = G.AllSubgroups().ToGroupWrapper();
            var names = NamesTree.BuildName(sg);
            G.Name = names[0].Name;
            DisplayGroup.HeadGenerators(G);
            Console.WriteLine(sg.Infos);
            names.Println("Group names");
            Console.WriteLine();
            Console.WriteLine();
            break;
        }
    }
}

(Dictionary<int, int> cycle, int state) Search(int n, Perm M, Perm Mr)
{
    var N = new Dictionary<int, int>();
    var infos = SearchStep(n, M, Mr, N);
    return infos;
}

(int ord, int state) ValidCycle(int dim, Dictionary<int, int> cycle)
{
    if (cycle.Values.Distinct().Count() != cycle.Count)
        return (-1, -1);

    var list = new List<int>() { 1 };
    var cyc = new Dictionary<int, int>(cycle);
    while (cyc.Count != 0)
    {
        var c0 = cyc.Keys.First();
        var cs = c0;
        var k = 0;
        do
        {
            ++k;
            if (!cyc.ContainsKey(c0))
                return (Lcm(list.ToArray()), -1);

            var c1 = cyc[c0];
            cyc.Remove(c0);
            c0 = c1;
        } while (cs != c0);

        list.Add(k);
    }

    return (Lcm(list.ToArray()), dim - cycle.Count);
}

(Dictionary<int, int> cycle, int state) SearchStep(int ord, Perm M, Perm Mr, Dictionary<int, int> N)
{
    var dim = M.Sn.N;
    var xn = dim.Range().ToHashSet();
    var remKeys = xn.Except(N.Keys).ToHashSet();
    var remValues = xn.Except(N.Values).ToHashSet();
    var N0 = new Dictionary<int, int>(N);

    var infos1 = ValidCycle(dim, N0);
    if (infos1.ord == ord)
    {
        foreach (var i in remKeys)
            N0[i] = i;

        return (N0, 1);
    }
    else if (infos1.state == -1)
    {
        if (infos1.ord == -1 || ord < infos1.ord || ord % infos1.ord != 0)
            return (N0, -1);
    }
    else if (N0.Count == dim)
        return (N0, -1);

    var attempts = new HashSet<int>();
    var a0 = remKeys.Min();
    while (true)
    {
        var a1 = a0;
        var setValues = remValues.Except(attempts).ToHashSet();
        if (setValues.Count == 0)
            return (N0, 0);

        var b1 = setValues.Max();
        attempts.Add(b1);
        var N1 = new Dictionary<int, int>(N0);
        N1[a1] = b1;

        var brType = -1;
        while (setValues.Any())
        {
            setValues.Remove(a1);
            var b2 = M.Table[b1];
            var a2 = Mr.Table[a1];
            if (N1.ContainsValue(b2) && !N1.ContainsKey(a2))
            {
                brType = 0;
                break;
            }

            if (N1.ContainsKey(a2))
            {
                brType = 1;
                break;
            }

            N1[a2] = b2;
            a1 = a2;
            b1 = b2;
        }

        if (brType == 0)
            continue;

        var infos3 = ValidCycle(dim, N1);
        if (infos3.ord == ord)
            return (N1, 1);
    }
}

{
    FrobPerm(4, 2, 3);
    FrobPerm(7, 3, 2);
    // FrobPerm(5, 10, 4);
    FrobPerm(5, 4, 2);
    // FrobPerm(5, 4, 3);
    FrobPerm(11, 5, 3);
}
