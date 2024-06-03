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
using FastGoat.UserGroup.Floats;
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
using GroupRegX = System.Text.RegularExpressions;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

void TestOrth<T>(ConcreteGroup<T> g, bool indStab = true, params (int dim, int[] linIdx)[] infos)
    where T : struct, IElt<T>
{
    var ct = FG.CharacterTableEmpty(g);
    ct.DerivedSubGroupLift();
    if (indStab)
        ct.InductionFromStabilizers();

    Logger.SetOff();
    DisplayGroup.HeadOrders(g);
    ct.DisplayCells(tableOnly: true);
    if (ct.TableComplete)
        return;

    ct.SolveSumSquare();
    ct.SolveOrthogonality(infos);
    ct.DisplayCells();
    Console.Beep();
}

void LinChisMulGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    Console.WriteLine(g.ShortName);
    var ct = FG.CharacterTable(g);
    if (!ct.TableComplete)
        ct.SolveOrthogonality();

    var chG = Group.MulGroup($"Ch({g})", ct.DoneChis.Where(chi => chi.IsLinear).ToArray());
    var gtype = Group.AbelianGroupType(chG);
    chG.Name = gtype.Glue(" x ", "C{0}");
    DisplayGroup.HeadElements(chG);
    DisplayGroup.Generators(chG);
}

void runOrth()
{
    GlobalStopWatch.Restart();
    TestOrth(FG.Quaternion(16), indStab: false);
    TestOrth(FG.MetaCyclicSdp(4, 4, 3), indStab: false);
    TestOrth(FG.DiCyclic(5), indStab: false, infos: [(2, [0, 2]), (2, [0, 2])]);
    TestOrth(FG.DiCyclic(6), indStab: false, infos: [(2, [0, 2]), (2, [0, 2])]);

    TestOrth(FG.SL2p(3));
    TestOrth(FG.GL2p(3));
    TestOrth(Group.SemiDirectProd(FG.SL2p(3), FG.AbelianMat(2)), infos: (2, [0, 1, 2, 3, 4, 5]));
    TestOrth(Product.Generate(FG.SL2p(3), FG.AbelianMat(2)), infos: (2, [0, 1, 2, 3, 4, 5]));

    GlobalStopWatch.Show(); // Time:54.696s
}

void KerChisGroup<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
{
    var ct = FG.CharacterTable(g);
    if (!ct.TableComplete)
        ct.SolveOrthogonality((2, [0, 1, 2]));

    ct.DisplayCells();

    var listKers = ct.AllCharacters.Select((chi, k) => (chi, k))
        .Select(e => Group.Generate($"Ker(Ꭓ.{e.k + 1})", g, e.chi.Kernel().ToArray()))
        .ToList();

    var start = new HashSet<ConcreteGroup<T>>(listKers, new GroupSetEquality<T>());
    var allNormals = new HashSet<ConcreteGroup<T>>(listKers, new GroupSetEquality<T>());
    var sz = 0;
    while (sz != allNormals.Count)
    {
        sz = allNormals.Count;
        foreach (var n1 in allNormals.ToArray())
        {
            foreach (var n2 in start)
            {
                if (n2.SubSetOf(n1) || n1.SubSetOf(n2))
                    continue;

                var n3 = Group.Generate($"{n1} ∩ {n2}", g, n1.Intersect(n2).ToArray());
                allNormals.Add(n3);
            }
        }
    }

    var gAllSubgrs = g.AllSubgroups();
    Console.WriteLine(gAllSubgrs.Infos);
    allNormals.OrderBy(n => n.Count()).Println(e => e.ShortName, $"Total:{allNormals.Count}");
    Console.WriteLine();
    if (gAllSubgrs.Infos.AllNorms != allNormals.Count)
        throw new();
}

ConcreteGroup<KMatrix<Cnf>> GetGLnC(ConcreteGroup<Mat> mtGL)
{
    var n = mtGL.Neutral().GL.N;
    var p = mtGL.Neutral().GL.P;

    var Up = FG.UnInt(p);
    var e0 = Up.GetGenerators().First();
    var cnf = Cnf.Nth(p - 1);
    var vals = mtGL.SelectMany(mat => mat.Table).Distinct().ToHashSet();
    var iso = (p - 1).Range().Select(k => (k, e0.Pow(k).K)).Where(e => vals.Contains(e.K))
        .ToDictionary(e => e.K, e => cnf.Pow(e.k).Simplify());
    iso[0] = cnf.Zero;
    return Group.MulGroup(mtGL.Name,
        mtGL.GetGenerators().Select(mat => mat.Table.Select(z => iso[z]).ToKMatrix(n)).ToArray());
}

Polynomial<Cnf, Xi>[] Prod(KMatrix<Cnf> A, Polynomial<Cnf, Xi>[] X)
{
    if (A.N != X.Length)
        throw new();

    var zero = X[0].Zero;
    return A.Rows.Select(row => row.Zip(X).Aggregate(zero, (sum, e) => sum + e.First * e.Second)).ToArray();
}

void test0()
{
    var d8gl = FG.DihedralGL2p(4);
    var d8 = GetGLnC(d8gl);
    DisplayGroup.HeadOrdersGenerators(d8);

    var (x0, y0, x1, y1) = Ring.Polynomial(Cnf.CnfZero, MonomOrder.Lex, "x0", "y0", "x1", "y1").Deconstruct();
    var i = Cnf.I;
    var X = new[] { x0.One, x0.One * 2 };
    Console.WriteLine("X");
    Console.WriteLine($"[{X.Glue(", ", "{0,10}")}]");
    var e4 = d8.First(e => d8.ElementsOrders[e] == 4);
    foreach (var A in Group.MulGroup("C4", e4))
    {
        Console.WriteLine("A");
        Console.WriteLine(A);
        Console.WriteLine("A*X");
        Console.WriteLine($"[{Prod(A, X).Glue(", ", "{0,10}")}]");
        Console.WriteLine();
    }
}

KMatrix<Cnf> Rotate2D(int k)
{
    var c = Cnf.Nth(k);
    return new[] { c.Re, -c.Im, c.Im, c.Re }.ToKMatrix(2);
}

KMatrix<Cnf> Rotate3DX(int k)
{
    var c = Cnf.Nth(k);
    var (z, o) = (c.Zero, c.One);
    return new[] { o, z, z, z, c.Re, -c.Im, z, c.Im, c.Re }.ToKMatrix(3);
}

KMatrix<Cnf> Rotate3DY(int k)
{
    var c = Cnf.Nth(k);
    var (z, o) = (c.Zero, c.One);
    return new[] { c.Re, z, -c.Im, z, o, z, c.Im, z, c.Re }.ToKMatrix(3);
}

KMatrix<Cnf> Rotate3DZ(int k)
{
    var c = Cnf.Nth(k);
    var (z, o) = (c.Zero, c.One);
    return new[] { c.Re, -c.Im, z, c.Im, c.Re, z, z, z, o }.ToKMatrix(3);
}

KMatrix<Cplx> CplxMatrix(KMatrix<Cnf> A) => A.Select(e => new Cplx(e.ToComplex)).ToKMatrix(A.M);

void test1()
{
    var g0 = FG.DicyclicGL2p(3);
    var g1 = GetGLnC(g0);

    var (z, o) = (Cnf.CnfZero, Cnf.CnfOne);
    var a0 = Rotate2D(3);
    var a1 = new[] { z, o, -o, z }.ToKMatrix(2);
    var g2 = Group.MulGroup("H", a0, a1);
    DisplayGroup.HeadElements(g2);
    DisplayGroup.AreIsomorphics(g0, g2);
    Ring.DisplayMatrix(Ring.Matrix(2,
        g2.GetGenerators().SelectMany(e => e).Select(e => FG.PrettyPrintCnf(e).c).ToArray()));
}

KMatrix<Cnf> RMatrix(KMatrix<Cnf> A, bool ignoreRealMatrix = true)
{
    if (A.M != A.N)
        throw new();

    var B = A.Select(c => c.Re).ToKMatrix(A.M);
    var C = A.Select(c => c.Im).ToKMatrix(A.M);

    if (!ignoreRealMatrix && C.IsZero())
        return A;

    Cnf Mij(int i, int j)
    {
        if (i < A.M)
            return j < A.M ? B[i, j] : C[i, j - A.M];
        else
            return j < A.M ? -C[i - A.M, j] : B[i - A.M, j - A.M];
    }

    return (2 * A.M).Range().Grid2D().Select(e => Mij(e.t1, e.t2)).ToKMatrix(2 * A.M);
}

ConcreteGroup<KMatrix<Cnf>> GetGLnR(ConcreteGroup<KMatrix<Cnf>> G)
{
    var gens = G.GetGenerators().Select(e => RMatrix(e)).ToArray();
    foreach (var e in gens)
    {
        Console.WriteLine(CplxMatrix(e));
        Console.WriteLine();
    }

    return Group.MulGroup(G.Name, gens);
}

KMatrix<Cnf> CrossProd(KMatrix<Cnf> U, KMatrix<Cnf> V)
{
    if (U.Dim != (3, 1) || V.Dim != (3, 1))
        throw new();

    var (u1, u2, u3) = U.Deconstruct();
    var (v1, v2, v3) = V.Deconstruct();
    return new[] { u2 * v3 - u3 * v2, u3 * v1 - u1 * v3, u1 * v2 - u2 * v1 }.ToKMatrix(3);
}

KPoly<Cnf> CharacPol(KMatrix<Cnf> M)
{
    var X = Ring.Polynomial(Cnf.CnfOne);
    var o = X.One;
    var M0 = M.Select(c => c * o).ToKMatrix(M.M);
    var I = M0.One;
    Console.WriteLine(M0 - X * I);
    var pol = Ring.Determinant((M0 - X * I).Coefs, o);
    return pol.ToKPoly(X.ExtractIndeterminate);
}

bool IsDiagPerm2D(KMatrix<Cnf> M)
{
    if (M.Dim != (2, 2))
        return false;

    var (e00, e01, e10, e11) = M.Select(c => c.IsZero()).Deconstruct();
    return (e00 && e11 && !e01 && !e10) || (!e00 && !e11 && e01 && e10);
}

KMatrix<Cnf>[] Change(KMatrix<Cnf>[] gens, int lcm0 = 0)
{
    var ns = gens.SelectMany(gen => gen).Select(c => c.N).Distinct().ToArray();
    var lcm = lcm0 == 0 ? Lcm(ns) : lcm0;
    // if (lcm % 4 != 0)
    //     lcm *= lcm % 2 == 0 ? 2 : 4;

    Console.WriteLine(lcm);
    return gens.Select(m => m.Select(c => c.ToCnfN(lcm)).ToKMatrix(m.M)).ToArray();
}

void test3()
{
    var g0 = FG.MetaCyclicSdpMat(7, 3, 2);
    var g1 = GetGLnC(g0);

    var (a0, a1) = g1.GetGenerators().Deconstruct();

    var P0 = CharacPol(a0).ToCnfPoly(21);
    var roots0 = IntFactorisation.AlgebraicRoots(P0.ToEPolyX());
    roots0.Println($"P0:{P0} F:{roots0[0].F}");

    var P1 = CharacPol(a1);
    var roots1 = IntFactorisation.AlgebraicRoots(P1.ToEPolyX());
    roots1.Println($"P:{P1} F:{roots1[0].F}");
}

