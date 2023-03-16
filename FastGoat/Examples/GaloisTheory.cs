using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class GaloisTheory
{
    public static ConcreteGroup<Perm> GaloisGroup<K>(List<EPoly<K>> roots) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (roots.ToHashSet().Count != roots.Count || roots.Select(e => e.F).Distinct().Count() != 1)
            throw new("Roots must all be differents and belong to the same extensions");

        var (X, a) = FG.EPolyXc(roots[0].F, 'a');
        var prod = roots.Select(r => X - r).Aggregate(X.One, (acc, xa) => acc * xa);
        var minPoly = a.F.Substitute(X);
        if (roots.Count != a.F.Degree || !roots.Select(r => X - r).Aggregate(X.One, (acc, xa) => acc * xa).Equals(a.F.Substitute(X)))
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

        Console.WriteLine($"Polynomial P = {minPoly}");
        roots.Select(r => X - r).Println("Factorization");
        Console.WriteLine();
        var Gal = Group.Generate("Gal( Q(a)/Q )", sn, sigmas.ToArray());
        DisplayGroup.HeadElements(Gal);

        return Gal;
    }

    public static void NormalExtensionCase()
    {
        var x = FG.QPoly();
        
        {
            var (X0, y0) = FG.EPolyXc(x.Pow(2) - 2, 'a');
            var gal = GaloisGroup(new List<EPoly<Rational>>() { y0, -y0 });
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2));
            Console.WriteLine();
        }

        {
            var (X0, y0) = FG.EPolyXc(x.Pow(2) + 4 * x - 2, 'a');
            var gal = GaloisGroup(new List<EPoly<Rational>>() { y0, -y0 - 4 });
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(3) - 3 * x - 1);
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(3));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) - 4 * x.Pow(2) + 2);
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(4));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) + 1);
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 2));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) - 4 * x.Pow(2) + 1);
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(2, 2));
            Console.WriteLine();
        }

        {
            var roots = IntFactorisation.SplittingField(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(4));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots);
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
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Symmetric(3));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + 243, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Symmetric(3));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + 12, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Symmetric(3));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + x.Pow(5) - 7 * x.Pow(4) - 2 * x.Pow(3) + 7 * x.Pow(2) + 2 * x - 1, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(6));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(6) + x.Pow(5) - 7 * x.Pow(4) - 2 * x.Pow(3) + 7 * x.Pow(2) + 2 * x - 1, 'y');
            var roots = IntFactorisation.AlgebraicFactors(y.F.Substitute(X)).Select(f => -f[0] / f[1]).ToList();
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(6));
            Console.WriteLine();
        }

        {
            var (X, y) = FG.EPolyXc(x.Pow(7) + x.Pow(6) - 18 * x.Pow(5) - 35 * x.Pow(4) + 38 * x.Pow(3) + 104 * x.Pow(2) + 7 * x - 49, 'y');
            var roots = IntFactorisation.AlgebraicRoots(y.F.Substitute(X));
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Abelian(6));
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

            var gal = GaloisGroup(roots);
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

            var gal = GaloisGroup(roots);
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

            var gal = GaloisGroup(roots);
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

            var gal = GaloisGroup(roots);
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

            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Quaternion(8));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) + 24 * x.Pow(4) + 16, 'y');
            var roots = IntFactorisation.AlgebraicRoots(r0.F.Substitute(X0));
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Dihedral(4));
            Console.WriteLine();
        }

        {
            var (X0, r0) = FG.EPolyXc(x.Pow(8) - 12 * x.Pow(6) + 36 * x.Pow(4) - 36 * x.Pow(2) + 9, 'y');
            var roots = IntFactorisation.AlgebraicRoots(r0.F.Substitute(X0), true);
            var gal = GaloisGroup(roots);
            DisplayGroup.AreIsomorphics(gal, FG.Quaternion(8));
            Console.WriteLine();
        }
    }
}