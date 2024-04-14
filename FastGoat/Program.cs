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

WordGroup[] GetGroups(int ord)
{
    return GroupExt.DB.Select(s => s.Split(';'))
        .Where(s => int.Parse(s[0]) == ord)
        .Select(s =>
        {
            Logger.Level = LogLevel.Off;
            var g = FG.WordGroup(s[1], s[2]);
            Logger.Level = LogLevel.Level1;
            return g;
        })
        .ToArray();
}

void MatrixFormGroupsOrd8()
{
    var gl2 = FG.GLnq(2, 5);
    var so3 = FG.GO3p(3);
    foreach (var g8 in GetGroups(8))
    {
        Console.WriteLine(g8.ShortName);
        if (g8.Name == "Q8")
        {
            var iso = FG.Quaternion(8);
            DisplayGroup.Generators(iso);
            DisplayGroup.Generators(Group.IsomorphicSubgroup(gl2, g8));
        }
        else
        {
            if (g8.Name != "C2 x C2 x C2")
                DisplayGroup.Generators(Group.IsomorphicSubgroup(gl2, g8));

            if (g8.Name != "C8")
                DisplayGroup.Generators(Group.IsomorphicSubgroup(so3, g8));
        }
    }
}

void MatrixFormGroupsOrd9()
{
    var gl2 = new GL(2, 19);
    var a = FG.FqX(19);
    var e3 = a.Pow(6)[0].K;
    var e9 = a.Pow(2)[0].K;
    var gen9 = gl2[e9, 0, 0, 1];
    var gen3a = gl2[e3, 0, 0, 1];
    var gen3b = gl2[1, 0, 0, e3];
    var g9 = Group.Generate("C9mat", gl2, gen9);
    var g33 = Group.Generate("(C3 x C3)mat", gl2, gen3a, gen3b);
    var ab9 = FG.Abelian(9);
    var ab33 = FG.Abelian(3, 3);

    DisplayGroup.Generators(g9);
    DisplayGroup.AreIsomorphics(g9, ab9);
    Console.WriteLine();

    DisplayGroup.Generators(g33);
    DisplayGroup.AreIsomorphics(g33, ab33);
    Console.WriteLine();
}

void MatrixFormGroupsOrd10()
{
    var gl2 = FG.GLnp(2, 11);
    foreach (var g10 in GetGroups(10))
    {
        var iso = Group.IsomorphicSubgroup(gl2, g10);
        DisplayGroup.Generators(iso);
    }
}

void MatrixFormGroupsOrd12()
{
    var gl2 = FG.GLnp(2, 7);
    DisplayGroup.HeadOrders(gl2);
    foreach (var g12 in GetGroups(12))
    {
        Console.WriteLine(g12.ShortName);
        if (g12.Name == "A4")
        {
            var so33 = FG.SO3q(3);
            var iso = Group.IsomorphicSubgroup(so33, g12);
            DisplayGroup.Generators(iso);
        }
        else
        {
            var iso = Group.IsomorphicSubgroup(gl2, g12);
            DisplayGroup.Generators(iso);
        }
    }
}

void MatrixFormGroupsOrd14()
{
    var gl2 = FG.GLnq(2, 8);
    foreach (var g14 in GetGroups(14))
    {
        var iso = Group.IsomorphicSubgroup(gl2, g14);
        DisplayGroup.Generators(iso);
    }
}

// {
//     MatrixFormGroupsOrd8();
//     MatrixFormGroupsOrd9();
//     MatrixFormGroupsOrd10();
//     MatrixFormGroupsOrd12();
//     MatrixFormGroupsOrd14();
// }

int OrderMatOrth(ZnInt x0, ZnInt y0)
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

(ZnInt x, ZnInt y)[] CandidateGO2p(ConcreteGroup<ZnInt> Zp, int n)
{
    var a = Zp.ElementsOrders.MaxBy(e => e.Value).Key;
    var square = Zp.Append(a.Zero).Select(x => (x, x2: x * x)).GroupBy(e => e.x2)
        .ToDictionary(e => e.Key, e => e.Select(f => f.x).ToArray());
    var dicSquare = Zp.ToDictionary(x => x, x => square.ContainsKey(x) ? square[x] : []);
    dicSquare[a.Zero] = [];
    return Zp.Append(a.Zero)
        .Select(x0 => (x: x0, yList: dicSquare[1 - x0 * x0]))
        .Where(e => e.yList.Length != 0)
        .SelectMany(e => e.yList.Select(y0 => (e.x, y: y0)))
        .Distinct()
        .Select(e => (e.x, e.y))
        .Where(e => OrderMatOrth(e.x, e.y) == n)
        .OrderBy(e => e.x.K)
        .ToArray();
}

(ZnInt x, ZnInt y) FindGOnp(int n)
{
    foreach (var p in Primes10000.Where(p => p >= n))
    {
        var a = Un.FirstGen(p);
        var Zp = Group.MulGroup($"F{p}", a);
        var candidates = CandidateGO2p(Zp, n);

        if (candidates.Length != 0)
            return candidates[0];
    }

    throw new();
}

(MatFq U, MatFq J) EigenVals(MatFq m)
{
    var z0 = m.GLnq.Fq.X.Zero;
    var (x, x0, x1) = Ring.Polynomial(z0, MonomOrder.Graded, "x", "x0", "x1").Deconstruct();
    var M = m.Table.Select(k => k * x.One).ToKMatrix(m.GLnq.N);
    Console.WriteLine(M);
    var M0 = (M - x * M.One);
    var P0 = M0[0, 0] * M0[1, 1] - M0[0, 1] * M0[1, 0];
    Console.WriteLine(P0);
    var P1 = P0.ToKPoly(x);

    var sols = IntFactorisation.FirrFsep(P1, m.GLnq.Fq.X);
    sols.Println($"P1 = {P1}");
    if (sols.Count != 2)
        throw new();
    
    var r0 = -sols[0].g[0] / sols[0].g[1];
    var r1 = -sols[1].g[0] / sols[1].g[1];

    var E = new[] { x0, x1 }.ToKMatrix(m.GLnq.N);
    KMatrix<EPoly<ZnInt>> E0, E1;

    var sys0 = M * E - r0 * x.One * E;
    Console.WriteLine(sys0);
    var sol0 = Ring.ReducedGrobnerBasis(sys0.Where(o => !o.IsZero()).ToArray());
    sol0.Println("Sys0");
    var (c00, c01) = (sol0[0][x0.ExtractMonom], sol0[0][x1.ExtractMonom]);
    if (!c00.IsZero())
    {
        if (!c01.IsZero())
            E0 = new[] { z0.One, (-c00 / c01) * z0.One }.ToKMatrix(m.GLnq.N);
        else
            E0 = new[] { z0, z0.One }.ToKMatrix(m.GLnq.N);
    }
    else
    {
        if (!c01.IsZero())
            E0 = new[] { z0.One, z0 }.ToKMatrix(m.GLnq.N);
        else
            throw new();
    }

    Console.WriteLine();

    var sys1 = M * E - r1 * x.One * E;
    Console.WriteLine(sys1);
    var sol1 = Ring.ReducedGrobnerBasis(sys1.Where(o => !o.IsZero()).ToArray());
    sol1.Println("Sys1");

    var (c10, c11) = (sol1[0][x0.ExtractMonom], sol1[0][x1.ExtractMonom]);
    if (!c10.IsZero())
    {
        if (!c11.IsZero())
            E1 = new[] { z0.One, (-c10 / c11) * z0.One }.ToKMatrix(m.GLnq.N);
        else
            E1 = new[] { z0, z0.One }.ToKMatrix(m.GLnq.N);
    }
    else
    {
        if (!c11.IsZero())
            E1 = new[] { z0.One, z0 }.ToKMatrix(m.GLnq.N);
        else
            throw new();
    }

    Console.WriteLine();

    var U = new[] { r0, z0, z0, r1 }.ToKMatrix(m.GLnq.N);
    var J = KMatrix<EPoly<ZnInt>>.MergeSameRows(E0, E1);
    Console.WriteLine("U");
    Console.WriteLine(U);
    Console.WriteLine("J");
    Console.WriteLine(J);
    Console.WriteLine("J * U * J-1");
    Console.WriteLine(J * U * J.Inv());
    Console.WriteLine();
    var gl = m.GLnq;
    return (U: gl.Create(U.Select(c => c).ToArray()), J: gl.Create(J.Select(c => c).ToArray()));
}

(int n, int p, ZnInt ord_2, ZnInt ord_n) FindP(int n)
{
    foreach (var p in Primes10000)
    {
        var a0 = Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
        if (a0 == -1)
            continue;

        var a = new ZnInt(p, a0);
        var Zp = Group.MulGroup($"F{p}", a);

        var ords_2_n = Zp.ElementsOrders.Where(e => e.Value == n || e.Value == 2)
            .GroupBy(e => e.Value)
            .ToDictionary(e => e.Key, e => e.Select(f => f.Key).ToArray());

        if (ords_2_n.Count != 2)
            continue;

        var ord_2 = ords_2_n[2].First();
        var ord_n_pow = ords_2_n[n].ToDictionary(e => e, e => e.Pow(n / 2));
        var ord_n = ords_2_n[n].OrderByDescending(e => !ord_n_pow[e].Equals(ord_2)).First();
        return (n, p, ord_2, ord_n);
    }

    throw new();
}

void FindGLnp(int n, int p, ZnInt ord_2, ZnInt ord_n)
{
    var gl = new GL(2, p);
    
    var a0 = Solve_k_pow_m_equal_one_mod_n_strict(p, p - 1);
    var a = new ZnInt(p, a0);
    var Zp = Group.MulGroup($"F{p}", a);

    var id = gl.Neutral();

    foreach (var c0 in Zp)
    {
        var p0 = gl[1, a.K, 0, 1];
        var p1 = gl[1, 0, c0.K, 1];
        var m0 = gl.Op(gl.Invert(p0), gl.Op(gl[ord_n.K, 0, 0, ord_n.Inv().K], p0));
        var m1 = gl.Op(gl.Invert(p1), gl.Op(gl[1, 0, 0, ord_2.K], p1));
        if (!gl.Times(gl.Op(m0, m1), 2).Equals(id))
            continue;
        
        var D2n = Group.Generate($"D{2 * n}", gl, m0, m1);
        DisplayGroup.Generators(D2n);
        var subgs = D2n.AllSubgroups();
        FG.FindIdGroup(D2n, subgs.Infos).Println(e => e.FullName);
        Console.WriteLine();
        return;
    }

    throw new();
}

void CaseD8()
{
    Ring.MatrixDisplayForm = Ring.MatrixDisplay.SquareBracketNoFmt;
    var d8 = Group.IsomorphicSubgroup(FG.GL2p(3), FG.Dihedral(4));
    DisplayGroup.HeadGenerators(d8);
    
    // Generators of D8
    // gen1 of order 2
    // [[1, 0],[2, 2]]
    // 
    // gen2 of order 4
    // [[2, 1],[1, 1]]
    // 

    var gl = new GLnq(2, 9);
    var a = gl[2, 1, 1, 1];
    var b = gl[1, 0, 2, 2];
    
    EigenVals(a);
    EigenVals(b);
}

{
    for (int n = 3; n < 33; ++n)
    {
        var (_, p, ord2, ordn) = FindP(n);
        FindGLnp(n, p, ord2, ordn);
    }
}