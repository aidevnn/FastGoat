using System.Numerics;
using System.Runtime.InteropServices.JavaScript;
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

EPoly<K>[] GetIndependants<K>(EPoly<K>[] list) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var bs = new List<EPoly<K>>();
    foreach (var a0 in list.Order())
    {
        if (bs.Count(a => a.Degree > 0) == 0)
        {
            if (a0.Degree > 0 || bs.Count == 0)
                bs.Add(a0);

            continue;
        }

        var q = bs.Where(a1 => !a1.IsZero())
            .Select(a1 => KMatrix<K>.MergeSameRows(new[] { a0, a1 }.Select(e => e.Poly.ToVMatrix(a0.F.Degree - 1)).ToArray()));
        if (q.All(mat => mat.NullSpace().nullity == 0))
            bs.Add(a0);
    }

    return bs.ToArray();
}

EPoly<K>[] GetBase<K>(EPoly<K> a) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    var bs = Array.Empty<EPoly<K>>();
    var n = a.F.Degree - 1;
    for (int i = 0; i <= n; ++i)
    {
        var ai = a.Pow(i);
        var bs0 = bs.Append(ai).ToArray();
        var mat = KMatrix<K>.MergeSameRows(bs0.Select(e => e.Poly.ToVMatrix(a.F.Degree - 1)).ToArray());
        if (mat.NullSpace().nullity == 0)
            bs = bs0.ToArray();
    }

    return bs;
}

void ExtDegree(EPoly<Rational> a)
{
    var bs = GetBase(a);
    Console.WriteLine($"[Q(a)/Q] = {bs.Length} with a={a}");
}

(EPoly<K> W, int l) Primitive3<K>(EPoly<K> U, EPoly<K> V) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    if (U.Degree == 0)
        return (V, 0);

    if (V.Degree == 0)
        return (U, 0);

    if (U.Degree < V.Degree)
        return Primitive3(V, U);

    var n = U.F.Degree;
    var vecV = V.Poly.ToVMatrix(n);

    for (int l = 1; l < 50; l++)
    {
        var W = -U - l * V;
        var M = KMatrix<K>.MergeSameRows(GetBase(W).Select(w => w.Poly.ToVMatrix(n)).ToArray());
        var vM = KMatrix<K>.MergeSameRows(vecV, M);
        var dimKerM = M.T.NullSpace().nullity;
        var dimKerVM = vM.T.NullSpace().nullity;
        if (dimKerM == dimKerVM)
        {
            return (W, l);
        }
    }

    throw new($"{U} ... {V}");
}

(KPoly<Rational> minPoly, EPoly<Rational> primElt) InvariantsSubGr(EPoly<Rational>[] roots, bool details = false)
{
    if (roots.Grid2D(roots).Any(e => !roots.Contains(e.t1.Substitute(e.t2))))
        throw new();

    var a = roots[0].X;
    var n = a.F.Degree;

    var (X, xis) = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, (n, "c"), "X");
    var R = xis.Select((c, i) => c * X.Pow(i)).Aggregate(X.Zero, (sum, xi) => sum + xi);
    var eqs = new List<Polynomial<Rational, Xi>>();
    foreach (var h in roots)
    {
        var Rh = xis.Select((c, i) => c * h.Pow(i).Poly.Substitute(X)).Aggregate(X.Zero, (sum, xi) => sum + xi);
        var e = Ring.Decompose(Rh - R, X.ExtractIndeterminate);
        eqs.AddRange(e.Item1.Values);
    }

    var mat = new KMatrix<Rational>(Rational.KZero(), eqs.Count, n);
    for (int i = 0; i < eqs.Count; i++)
    {
        var eq = eqs[i];
        for (int j = 0; j < n; j++)
        {
            mat.Coefs[i, j] = eq[xis[j].ExtractMonom];
        }
    }

    var sols = mat.NullSpace().Item2;
    
    var Rs = sols.Cols.Select(m => n.Range().Aggregate(a.Zero, (sum, i) => m[i, 0] * a.Pow(i) + sum))
        .Where(p => !p.IsZero()).Order().ToArray();

    // var primElt = Rs.Where(s0 => s0.Degree != 0).Aggregate(a.Zero, (e, s1) => Primitive3(e, s1).W);
    var primElt = Rs.Length == 1 ? Rs[0] : Rs.Where(p => p.Degree > 0).Min();

    if(details)
    {
        Console.WriteLine("System");
        Console.WriteLine(mat);
        Console.WriteLine("Sols");
        Console.WriteLine(sols);
        Console.WriteLine();
        Rs.Println($"Invariants");
    }

    foreach (var b in roots)
    {
        var coefs = n.Range().Select(i => (i, c: Rng.Next(-9, 10) * Rational.KOne())).ToArray();
        var re = coefs.Aggregate(primElt.Zero, (acc, e) => primElt.Pow(e.i) * e.c + acc); // random element x
        var f_re = coefs.Aggregate(primElt.Zero, (acc, e) => primElt.Substitute(b).Pow(e.i) * e.c + acc); // F(x)

        var str_x = new KPoly<Rational>('q', Rational.KOne(), coefs.Select(e => e.c).ToArray()).ToString();
        var str_bx = str_x.Replace("q", "F(q)");

        if (details)
        {
            Rs.Select(e => (e, e.Substitute(b).Equals(e))).Println();
            
            Console.WriteLine($"PrimElt q={primElt} ");
            Console.WriteLine($"Is [F(q) = q] {primElt.Substitute(b).Equals(primElt)}");
            Console.WriteLine("Random Test");
            Console.WriteLine($"  x  = s(q)    = {str_x}");
            Console.WriteLine($"F(x) = s(F(q)) = {str_bx}");
            Console.WriteLine($"  x  = s(q)    = {re}");
            Console.WriteLine($"F(x) = F(s(q)) = {re.Substitute(b)}");
            Console.WriteLine($"F(x) = s(F(q)) = {f_re}");
            Console.WriteLine();
        }

        if (!re.Equals(re.Substitute(b)) || !re.Equals(f_re) || Rs.Any(e => !e.Substitute(b).Equals(e)))
            throw new();
    }

    var X0 = FG.KPoly('X', a);
    var r = IntFactorisation.Norm(X0 - primElt);
    var sff = IntFactorisation.YunSFF(r);
    if (sff.Count != 1 || sff[0].q != 1)
        throw new();

    var sep = sff[0];
    Console.WriteLine($"MinPoly = {sep}");
    ExtDegree(primElt);
    Console.WriteLine();
    return (sep.g, primElt);
}

HashSet<T> DepthFirstSearch<T>(T[] nodes, ILookup<T, T> edges) where T : IEquatable<T>
{
    var a0 = nodes.First();
    var visited = new HashSet<T>() { a0 };
    var stack = new Stack<T>(new[] { a0 });

    while (stack.Count > 0)
    {
        var ai = stack.Pop();
        var adjacencies = edges[ai];
        foreach (var aj in adjacencies.Except(visited))
        {
            stack.Push(aj);
            visited.Add(aj);
        }
    }

    return visited;
}

(KPoly<Rational> minPoly, EPoly<Rational> primElt) Inter(KPoly<Rational> Pa, KPoly<Rational> Pb)
{
    var (Xa, a0) = FG.EPolyXc(Pa, 'a');
    var (Xb, b0) = FG.EPolyXc(Pb, 'b');
    var P1 = IntFactorisation.AlgebraicFactors(Pb.Substitute(Xa), true)[0];
    var pea = IntFactorisation.PrimitiveElt(Pa, Pb)[0];

    var Pc = pea.F;
    var (Xc, c) = FG.EPolyXc(Pc, 'c');

    Console.WriteLine(Pc);
    var Pca = IntFactorisation.AlgebraicFactors(Pc.Substitute(Xa), true);
    var Pcb = IntFactorisation.AlgebraicFactors(Pc.Substitute(Xb), true);

    var ac = pea.a.Substitute(c);
    var bc = pea.b.Substitute(c);
    var l = (c - bc) / ac;
    var P1c = P1.SubstituteP0b(Xc, ac).Substitute(Xc - l * ac);

    Console.WriteLine();
    Console.WriteLine(l);
    Console.WriteLine(P1);
    Console.WriteLine(P1c);
    Console.WriteLine((Pa.Substitute(ac), Pb.Substitute(bc), P1c.Substitute(c)));

    var Ai = Pca.Select(p => p.SubstituteP0b(Xc, ac)).ToArray();
    var Bj = Pcb.Select(p => p.SubstituteP0b(Xc, bc)).ToArray();

    var nodes = Ai.Concat(Bj).ToArray();
    var edges = nodes.Grid2D(nodes).Where(e => !e.t1.Equals(e.t2) && Ring.Gcd(e.t1, e.t2).Degree > 0).ToLookup(e => e.t1, e => e.t2);
    nodes.Println();
    var dfs = DepthFirstSearch(nodes, edges);
    dfs.Println();
    var Px = dfs.Aggregate((e0, e1) => e0 * e1);
    Console.WriteLine(new { Px });
    var coefs = Px.Coefs;
    var qc = GetIndependants(coefs).Aggregate(c.Zero, (qci, ci) => Primitive3(qci, ci).W);
    GetIndependants(coefs).Println();
    Console.WriteLine(qc);
    ExtDegree(qc);

    var r = IntFactorisation.Norm(Xc - qc);
    var sff = IntFactorisation.YunSFF(r);
    if (sff.Count != 1 || sff[0].q != 1)
        throw new();

    var minPoly = sff[0].g;
    return (minPoly, qc);
}

(KPoly<Rational> minPoly, EPoly<Rational> primElt) Inter1(EPoly<Rational> a, EPoly<Rational> b)
{
    var Pa = IntFactorisation.CharacPoly(a, true);
    var Pb = IntFactorisation.CharacPoly(b, true);
    return Inter(Pa, Pb);
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    var x = FG.QPoly('X');

    var P1 = x.Pow(2) + 1;
    var P2 = x.Pow(4) + 3 * x.Pow(2) + 3;
    var P3 = x.Pow(6) + 12;
    var P4 = x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1;
    var P5 = x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1;
    var P6 = x.Pow(8) - 12 * x.Pow(6) + 23 * x.Pow(4) - 12 * x.Pow(2) + 1;
    var P7 = x.Pow(8) + 4 * x.Pow(6) + 2 * x.Pow(4) + 28 * x.Pow(2) + 1;
    var P8 = x.Pow(4) - x.Pow(2) + 1;
    var P9 = x.Pow(8) - x.Pow(4) + 1;
    var P10 = x.Pow(8) + 28 * x.Pow(4) + 2500;

    var roots = IntFactorisation.AlgebraicRoots(P6, true);
    var gal = GaloisTheory.GaloisGroup(roots);

    var allSubs = gal.Grid2D(gal).Select((e, k) => Group.Generate($"G{k}", gal, new[] { e.t1, e.t2 }))
        .OrderByDescending(g0 => g0.Count()).ToHashSet(new GroupSetEquality<Perm>());

    var sn = new Sn(roots.Count);
    var idxRoots = roots.Select((c0, k) => (k, c0)).ToDictionary(e => e.c0, e => e.k);
    var perm2roots = roots.Select(c0 => (c0, sn.CreateElement(roots.Select(c1 => idxRoots[c1.Substitute(c0)] + 1).ToArray())))
        .ToDictionary(e => e.Item2, e => e.c0);

    foreach (var sub in allSubs)
    {
        var rs = sub.Select(e => perm2roots[e]).ToArray();
        if (rs.Grid2D(rs).Any(e => !rs.Contains(e.t1.Substitute(e.t2))))
            throw new();

        Console.WriteLine("#####################################");
        InvariantsSubGr(rs, details: false);
        DisplayGroup.Head(sub);
    }
}