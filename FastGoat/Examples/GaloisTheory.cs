using System.Diagnostics.CodeAnalysis;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public record GaloisCorrespondence(
    ConcreteGroup<Perm> SubGr, 
    EPoly<Rational>[] roots, 
    EPoly<Rational> primElt,
    KPoly<Rational> minPoly);

public static class GaloisTheory
{
    public static ConcreteGroup<Perm> GaloisGroup<K>(List<EPoly<K>> roots, char primEltChar = 'α', bool details = false)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (roots.ToHashSet().Count != roots.Count || roots.Select(e => e.F).Distinct().Count() != 1)
            throw new("Roots must all be differents and belong to the same extensions");

        var (X, a) = FG.EPolyXc(roots[0].F, primEltChar);
        var prod = roots.Select(r => X - r).Aggregate(X.One, (acc, xa) => acc * xa);
        if (roots.Count != a.F.Degree || !prod.Equals(a.F.Substitute(X)))
            throw new("Extension must be normal");

        var n = a.F.Degree;
        var Fi = roots.Select((k, i) => (k, i)).ToDictionary(e => e.i, e => e.k.Poly);
        var idx = roots.Select((k, i) => (k, i)).ToDictionary(e => e.k, e => e.i);

        var sn = new Sn(n);
        var sigmas = new List<Perm>();
        for (int j = 0; j < n; j++)
        {
            var s_j = new List<int>();
            for (int i = 0; i < n; i++)
            {
                var asj = Fi[i].Substitute(roots[j]);
                s_j.Add(idx[asj]);
            }

            var as_j = sn.CreateElement(s_j.Select(k => k + 1).ToArray());
            sigmas.Add(as_j);
        }

        var Gal = Group.Generate($"Gal( Q({primEltChar})/Q )", sn, sigmas.ToArray());

        if (details)
        {
            var minPoly = a.F.Substitute(X);
            Console.WriteLine($"Polynomial P = {minPoly}");
            roots.Select(r => X - r).Println($"Factorization in Q({primEltChar})[X] with P({primEltChar}) = 0");
            Console.WriteLine();
            DisplayGroup.HeadElements(Gal);
        }

        return Gal;
    }

    public static EPoly<K>[] GetBase<K>(EPoly<K> a) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var bs = Array.Empty<EPoly<K>>();
        var n = a.F.Degree - 1;
        for (int i = 0; i <= n; ++i)
        {
            var ai = a.Pow(i);
            var bs0 = bs.Append(ai).ToArray();
            var mat = KMatrix<K>.MergeSameRows(bs0.Select(e => e.Poly.ToVMatrix(n)).ToArray());
            if (mat.NullSpace().nullity == 0)
                bs = bs0.ToArray();
        }

        return bs;
    }

    public static (EPoly<K> W, int l) PrimitiveEltComb<K>(EPoly<K> U, EPoly<K> V) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
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

        throw new($"U={U} and V={V}");
    }

    public static KPoly<K> Rewrite<K>(EPoly<K> a, EPoly<K> b) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (!a.F.Equals(b.F))
            throw new("Elements must belong to the same field");

        var n = b.F.Degree - 1;
        var bs = GetBase(a);
        var mat = KMatrix<K>.MergeSameRows(bs.Append(b).Select(e => e.Poly.ToVMatrix(n)).ToArray());
        var ns = mat.NullSpace();
        if (ns.nullity != 0)
        {
            return Ring.ReducedRowsEchelonForm(mat).A0.Cols.Last().ToKPoly('a');
        }

        return a.Zero.Poly.SubstituteChar('a');
    }

    public static void ExtDegree(EPoly<Rational> a)
    {
        var bs = GetBase(a);
        Console.WriteLine($"[Q(a)/Q] = {bs.Length} with a={a}");
    }

    public static (KPoly<Rational> minPoly, EPoly<Rational> primElt) InvariantsField(EPoly<Rational>[] roots, bool details = false)
    {
        if (roots.Grid2D(roots).Any(e => !roots.Contains(e.t1.Substitute(e.t2))))
            throw new();

        var a = roots[0].X;
        var n = a.F.Degree;

        var rows = new HashSet<KMatrix<Rational>>();
        foreach (var h in roots)
        {
            var hi = n.Range().Select(i => h.Pow(i)).ToArray();
            var h0i = n.Range().Select(i => hi.Select((h0, j) => h0[i] - (i == j ? 1 : 0)).ToKMatrix()).ToArray();
            rows.UnionWith(h0i);
        }

        var mat = KMatrix<Rational>.MergeSameCols(rows.ToArray());
        var sols = mat.NullSpace().Item2;

        var Rs = sols.Cols.Select(m => n.Range().Aggregate(a.Zero, (sum, i) => m[i, 0] * a.Pow(i) + sum)).Order().ToArray();

        if (details)
        {
            Console.WriteLine("System");
            Console.WriteLine(mat);
            Console.WriteLine("Sols");
            Console.WriteLine(sols);
            Console.WriteLine();
            Rs.Println($"Invariants");
        }

        var primElt = Rs.First(e => GetBase(e).Length * roots.Length == n);
        foreach (var b in roots)
        {
            var coefs = n.Range().Select(i => (i, c: IntExt.Rng.Next(-9, 10) * Rational.KOne())).ToArray();
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

        if (details)
        {
            Console.WriteLine($"MinPoly = {sep}");
            ExtDegree(primElt);
            Console.WriteLine();
        }

        return (sep.g, primElt);
    }

    public static IEnumerable<GaloisCorrespondence> SubFields(KPoly<Rational> P, int nbGens = 2, char primEltChar = 'y',
        bool details = false)
    {
        var roots = IntFactorisation.AlgebraicRoots(P, details: details);
        return SubFields(roots, nbGens, primEltChar, details);
    }
    
    public static IEnumerable<GaloisCorrespondence> SubFields(List<EPoly<Rational>> roots, int nbGens = 2, char primEltChar = 'y', bool details = false)
    {
        var P = roots[0].F;
        var gal = GaloisGroup(roots, primEltChar, details);

        var allSubs = nbGens == 2
            ? gal.Grid2D(gal).Select((e, k) => Group.Generate($"G{k}", gal, new[] { e.t1, e.t2 }))
                .OrderByDescending(g0 => g0.Count()).ToHashSet(new GroupSetEquality<Perm>())
            : gal.Grid3D(gal, gal).Select((e, k) => Group.Generate($"G{k}", gal, new[] { e.t1, e.t2, e.t3 }))
                .OrderByDescending(g0 => g0.Count()).ToHashSet(new GroupSetEquality<Perm>());

        var sn = new Sn(roots.Count);
        var idxRoots = roots.Select((c0, k) => (k, c0)).ToDictionary(e => e.c0, e => e.k);
        var perm2roots = roots.Select(c0 => (c0, sn.CreateElement(roots.Select(c1 => idxRoots[c1.Substitute(c0)] + 1).ToArray())))
            .ToDictionary(e => e.Item2, e => e.c0);

        var i = 1;
        foreach (var subGr in allSubs)
        {
            subGr.SetName($"G{i++}");

            var rs = subGr.Select(e => perm2roots[e]).ToArray();
            if (rs.Grid2D(rs).Any(e => !rs.Contains(e.t1.Substitute(e.t2))))
                throw new();

            if (details)
                Console.WriteLine("#####################################");
            var inv = GaloisTheory.InvariantsField(rs, details);
            if (details)
                DisplayGroup.Head(subGr);

            yield return new(subGr, rs, inv.primElt, inv.minPoly);
        }

        if (details)
        {
            Console.WriteLine($"############# End Galois Correspondence for polynomial P={P}");
            Console.WriteLine();
        }
    }

    public static void NormalExtensionCase()
    {
        var x = FG.QPoly();

        {
            var (X0, y0) = FG.EPolyXc(x.Pow(2) - 2, 'a');
            var gal = GaloisGroup(new List<EPoly<Rational>>() { y0, -y0 }, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2));
            Console.WriteLine();
        }

        {
            var (X0, y0) = FG.EPolyXc(x.Pow(2) + 4 * x - 2, 'a');
            var gal = GaloisGroup(new List<EPoly<Rational>>() { y0, -y0 - 4 }, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(3) - 3 * x - 1);
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(3));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) - 4 * x.Pow(2) + 2);
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(4));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) + 1);
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 2));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) - 4 * x.Pow(2) + 1);
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 2));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(4));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(5));
            Console.WriteLine();
        }

        {
            var (X0, y0) = FG.EPolyXc(x.Pow(6) + 243, 'a');

            var z = y0.Pow(3) / 9;
            var r0 = (-1 + z) / 2;
            var r1 = (-1 - z) / 2;
            var r2 = (1 + z) / 2;
            var r3 = (1 - z) / 2;

            var roots = new List<EPoly<Rational>>() { y0, -y0, r0 * y0, r1 * y0, r2 * y0, r3 * y0 };
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Symmetric(3));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + 243, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Symmetric(3));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + 12, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Symmetric(3));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + x.Pow(5) - 7 * x.Pow(4) - 2 * x.Pow(3) + 7 * x.Pow(2) + 2 * x - 1, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(6));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + x.Pow(5) - 7 * x.Pow(4) - 2 * x.Pow(3) + 7 * x.Pow(2) + 2 * x - 1, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(6));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(7) + x.Pow(6) - 18 * x.Pow(5) - 35 * x.Pow(4) + 38 * x.Pow(3) + 104 * x.Pow(2) + 7 * x - 49,
                'y');
            var roots = IntFactorisation.AlgebraicRoots(y.F.Substitute(X));
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(7));
            Console.WriteLine();
        }
    }

    public static void NormalExtensionCaseOrder8()
    {
        var x = FG.QPoly();

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) - 8 * x.Pow(6) + 20 * x.Pow(4) - 16 * x.Pow(2) + 2, 'y');

            var roots2 = IntFactorisation.AlgebraicFactors(X0.Pow(4) - 8 * X0.Pow(3) + 20 * X0.Pow(2) - 16 * X0 + 2)
                .Select(f => -f[0] / f[1]).ToList();

            var roots = new List<EPoly<Rational>>();
            foreach (var y in roots2)
            {
                roots.AddRange(IntFactorisation.AlgebraicFactors(X0.Pow(2) - y).Select(f => -f[0] / f[1]));
            }

            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(8));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) + 1, 'y');

            var roots2 = IntFactorisation.AlgebraicFactors(X0.Pow(4) + 1).Select(f => -f[0] / f[1]).ToList();

            var roots = new List<EPoly<Rational>>();
            foreach (var y in roots2)
            {
                roots.AddRange(IntFactorisation.AlgebraicFactors(X0.Pow(2) - y).Select(f => -f[0] / f[1]));
            }

            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 4));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) - x.Pow(4) + 1, 'y');

            var roots2 = IntFactorisation.AlgebraicFactors(X0.Pow(4) - X0.Pow(2) + 1).Select(f => -f[0] / f[1]).ToList();

            var roots = new List<EPoly<Rational>>();
            foreach (var y in roots2)
            {
                roots.AddRange(IntFactorisation.AlgebraicFactors(X0.Pow(2) - y).Select(f => -f[0] / f[1]));
            }

            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 2, 2));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) + 24 * x.Pow(4) + 16, 'y');

            var roots2 = IntFactorisation.AlgebraicFactors(X0.Pow(4) + 24 * X0.Pow(2) + 16).Select(f => -f[0] / f[1]).ToList();

            var roots = new List<EPoly<Rational>>();
            foreach (var y in roots2)
            {
                roots.AddRange(IntFactorisation.AlgebraicFactors(X0.Pow(2) - y).Select(f => -f[0] / f[1]));
            }

            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Dihedral(4));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9, 'y');

            var roots2 = IntFactorisation.AlgebraicFactors(X0.Pow(4) - 12 * X0.Pow(3) + 36 * X0.Pow(2) - 36 * X0 + 9)
                .Select(f => -f[0] / f[1]).ToList();

            var roots = new List<EPoly<Rational>>();
            foreach (var y in roots2)
            {
                roots.AddRange(IntFactorisation.AlgebraicFactors(X0.Pow(2) - y).Select(f => -f[0] / f[1]));
            }

            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Quaternion(8));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) + 24 * x.Pow(4) + 16, 'y');
            var roots = IntFactorisation.AlgebraicRoots(r0.F.Substitute(X0));
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Dihedral(4));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9, 'y');
            var roots = IntFactorisation.AlgebraicRoots(r0.F.Substitute(X0));
            var gal = GaloisGroup(roots, details: true);
            DisplayGroup.AreIsomorphics(gal, FG.Quaternion(8));
            Console.WriteLine();
        }
    }

    static void GaloisCorrespondence(KPoly<Rational> P, int nbGens = 2)
    {
        var subFields = SubFields(P, nbGens, details: true).ToArray();
    }

    public static void GaloisCorrespondenceSubfieldsExamples()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var x = FG.QPoly('X');

        GaloisCorrespondence(x.Pow(2) + 1);
        GaloisCorrespondence(x.Pow(4) - 10 * x.Pow(2) + 1);
        GaloisCorrespondence(x.Pow(4) - x.Pow(2) + 1);
        GaloisCorrespondence(x.Pow(6) + 12);
        GaloisCorrespondence(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
        GaloisCorrespondence(x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1); // Simple PGroup, no subField
        
        GaloisCorrespondence(x.Pow(8) - 12 * x.Pow(6) + 23 * x.Pow(4) - 12 * x.Pow(2) + 1, nbGens:3);
        GaloisCorrespondence(x.Pow(8) + 4 * x.Pow(6) + 2 * x.Pow(4) + 28 * x.Pow(2) + 1);
        GaloisCorrespondence(x.Pow(8) - x.Pow(4) + 1, nbGens:3);
        GaloisCorrespondence(x.Pow(8) + 28 * x.Pow(4) + 2500);
        
        GaloisCorrespondence(x.Pow(6) - 30 * x.Pow(4) + 225 * x.Pow(2) + 823);
        GaloisCorrespondence(x.Pow(6) + x.Pow(5) - 7 * x.Pow(4) - 2 * x.Pow(3) + 7 * x.Pow(2) + 2 * x - 1);
        
        GaloisCorrespondence(x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9);
        // GaloisCorrespondence(x.Pow(10) - 2 * x.Pow(9) - 20 * x.Pow(8) + 2 * x.Pow(7) + 69 * x.Pow(6) - x.Pow(5) - 69 * x.Pow(4) +
        //     2 * x.Pow(3) + 20 * x.Pow(2) - 2 * x - 1); // Dihedral 10
    }
}