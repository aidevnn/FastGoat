using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");

(EPoly<K> W, int l) Primitive3<K>(EPoly<K> U, EPoly<K> V) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    if (U.Degree == 0)
        return (V, 0);

    if (V.Degree == 0)
        return (U, 0);

    var n = U.F.Degree;
    var vecV = V.Poly.ToVMatrix(n);

    for (int l = 1; l < 50; l++)
    {
        var W = U + l * V;
        var M = KMatrix<K>.MergeSameRows(n.Range().Select(i => W.Pow(i).Poly.ToVMatrix(n)).ToArray());
        var vM = KMatrix<K>.MergeSameRows(vecV, M);
        var dimKerM = M.T.NullSpace().nullity;
        var dimKerVM = vM.T.NullSpace().nullity;
        if (dimKerM == dimKerVM)
        {
            return (W, l);
        }
    }

    throw new($"{U} / {V}");
}

void Invariants(EPoly<Rational> b)
{
    var a = b.X;
    var n = a.F.Degree;

    var eqs = KMatrix<Rational>.MergeSameRows(n.Range().Select(k => (b.Pow(k) - a.Pow(k)).Poly.ToVMatrix(n - 1)).ToArray());

    var ns = eqs.NullSpace();
    Console.WriteLine(ns.Item2);
    var Rs = ns.Item2.Cols.Select(m => n.Range().Aggregate(a.Zero, (sum, i) => m[i, 0] * a.Pow(i) + sum))
        .Where(p => !p.IsZero()).Order().ToArray();
    Rs.Println($"Invariants of {b.F.x} -> F({b.F.x}) = {b}");

    Rs.Select(e => (e, e.Substitute(b).Equals(e))).Println();
    var X = FG.KPoly('X', a);
    var q = Rs.Where(s0 => s0.Degree != 0).Aggregate(a.Zero, (e, s1) => Primitive3(e, s1).W);
    Console.WriteLine($"PrimElt q={q} ");
    Console.WriteLine($"Is [F(q) = q] {q.Substitute(b).Equals(q)}");
    IntFactorisation.NormDetails(X - q);
    Console.WriteLine();
}

{
    Ring.DisplayPolynomial = MonomDisplay.Caret;
    var x = FG.QPoly('X');

    var P1 = x.Pow(2) + 1;
    var P2 = x.Pow(4) + 3 * x.Pow(2) + 3;
    var P3 = x.Pow(6) + 12;
    var P4 = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
    var P5 = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
    var P6 = x.Pow(8) - 12 * x.Pow(6) + 23 * x.Pow(4) - 12 * x.Pow(2) + 1;
    var P7 = x.Pow(8) + 4 * x.Pow(6) + 2 * x.Pow(4) + 28 * x.Pow(2) + 1;;
    // var P = x.Pow(8) - x.Pow(4) + 1;
    // var P = x.Pow(8) + 28 * x.Pow(4) + 2500;
    // var P = x.Pow(4) - 2 * x.Pow(2) + 9;
    var roots = IntFactorisation.AlgebraicRoots(P7, details: true);
    var gal = GaloisTheory.GaloisGroup(roots);

    var n = roots.Count;
    var allSubs = gal.Grid2D(gal).Select((e, k) => Group.Generate($"G{k}", gal, new[] { e.t1, e.t2 }))
        .OrderByDescending(g0 => g0.Count()).ToHashSet(new GroupSetEquality<Perm>());

    var sn = new Sn(n);
    Perm.Style = DisplayPerm.Table;
    var idxRoots = roots.Select((c0, k) => (k, c0)).ToDictionary(e => e.c0, e => e.k);
    var perm2roots = roots.Select(c0 => (c0, sn.CreateElement(roots.Select(c1 => idxRoots[c1.Substitute(c0)] + 1).ToArray())))
        .ToDictionary(e => e.Item2, e => e.c0);

    foreach (var r0 in roots)
    {
        Invariants(r0);
    }
    
    foreach (var g0 in allSubs)
    {
        DisplayGroup.HeadElements(g0);
        var rs = g0.Select(e => perm2roots[e]).ToHashSet();
        if (rs.Grid2D(rs).Any(e => !rs.Contains(e.t1.Substitute(e.t2))))
            throw new();
    }
}
