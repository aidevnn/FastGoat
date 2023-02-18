using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class GaloisTheory
{
    public static void GaloisGroup<K>(List<EPoly<K>> roots) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
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
    }

    public static void NormalExtensionCase()
    {
        {
            var x = FG.QPoly();
            var (X0, y0) = FG.EPolyXc(x.Pow(2) - 2, 'a');
            GaloisGroup(new List<EPoly<Rational>>() { y0, -y0 });
        }

        {
            var x = FG.QPoly();
            var (X0, y0) = FG.EPolyXc(x.Pow(2) + 4 * x - 2, 'a');
            GaloisGroup(new List<EPoly<Rational>>() { y0, -y0 - 4 });
        }

        {
            var x = FG.QPoly();
            var roots = IntFactorisation.SplittingField(x.Pow(3) - 3 * x - 1);
            GaloisGroup(roots);
        }

        {
            var x = FG.QPoly();
            var roots = IntFactorisation.SplittingField(x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1);
            GaloisGroup(roots);
        }

        {
            var x = FG.QPoly();
            var roots = IntFactorisation.SplittingField(x.Pow(4) - 4 * x.Pow(2) + 2);
            GaloisGroup(roots);
        }

        {
            var x = FG.QPoly();
            var (X0, y0) = FG.EPolyXc(x.Pow(6) + 243, 'a');

            var z = y0.Pow(3) / 9;
            var r0 = (-1 + z) / 2;
            var r1 = (-1 - z) / 2;
            var r2 = (1 + z) / 2;
            var r3 = (1 - z) / 2;

            var roots = new List<EPoly<Rational>>() { y0, -y0, r0 * y0, r1 * y0, r2 * y0, r3 * y0 };
            GaloisGroup(roots);
        }
    }
}