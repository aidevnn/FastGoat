using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class PolynomialFactorization
{
    private static Random rnd = new Random();

    static KPoly<K> RandPoly<K>(K scalar, int p, int n) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var coefs = (n + 1).Range().Select(i => rnd.Next(-p, p + 1) * scalar.One).TrimSeq().ToArray();
        return new KPoly<K>('x', scalar, coefs);
    }

    static KPoly<ZnInt> ProdIrr(int p, int d)
    {
        var x = FG.ZPoly(p);
        return x.Pow(p.Pow(d)) - x;
    }

    static bool IsIrreductibleFp(KPoly<ZnInt> f)
    {
        if (f.Degree < 1)
            return true;

        var p = f.P;
        var n = f.Degree;
        var divs = IntExt.Dividors(n);
        return ProdIrr(p, n).Div(f).rem.IsZero() && divs.All(d => !ProdIrr(p, d).Div(f).rem.IsZero());
    }

    public static void IrreductibleRandPolys()
    {
        var n = 5;
        var p = 5;
        for (int i = 0; i < 4 * n; i++)
        {
            var f = RandPoly(ZnInt.KZero(p), p, n);
            Console.WriteLine($"{{0,{-7 * (n + 1)}}} is irreductible : {{1}}", f, IsIrreductibleFp(f));
        }
    }

    static List<(KPoly<K> g, int q, int i)> MusserSFF<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var L = new List<(KPoly<K>, int, int)>();
        var c = Ring.Gcd(f, f.Derivative);
        var i = 1;
        var g = f / c;
        while (g.Degree >= 1)
        {
            var p = Ring.Gcd(c, g);
            c = c / p;
            if (g.Degree > p.Degree)
                L.Add(((g / p).Monic, 1, i));

            g = p;
            ++i;
        }

        return L;
    }

    static List<(KPoly<K> g, int q, int i)> YunSFF<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var L = new List<(KPoly<K>, int, int)>();
        var l = 1;
        var df = f.Derivative;
        var u = Ring.Gcd(f, df);
        var v = f / u;
        var w = f.Derivative / u;
        while (v.Degree >= 1)
        {
            var w_dv = w - v.Derivative;
            var h = Ring.Gcd(v, w_dv);
            w = w_dv / h;
            v = v / h;
            if (h.Degree >= 1)
                L.Add((h.Monic, 1, l));

            ++l;
        }

        return L;
    }

    static (KPoly<K> c, int q) DeflateP<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var p = f.KZero.P;
        if (p == 0)
            throw new ArgumentException();

        var d = f.Degree;
        var decomp = IntExt.PrimesDecomposition(d).ToArray();
        if (!decomp.Contains(p))
            return (f, 0);

        var q = decomp.Count(i => i == p);
        var xi = q.Range(1).Select(i => p.Pow(i)).ToArray();
        var dq = (d + 1).Range();
        if (dq.Where(j => j % p != 0).All(j => f[j].IsZero()))
        {
            var coefs = dq.Where(j => j % p == 0).Select(j => f[j]).ToArray();
            return (new KPoly<K>(f.x, f.KZero, coefs), p);
        }

        return (f, 0);
    }

    static IEnumerable<(KPoly<K> g, int q, int m)> GianniTrager<K>(KPoly<K> f, int q)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (f.Degree != 0)
        {
            var sff = MusserSFF(f);
            var l0 = sff.Select(e => (e.g, q, m: e.i)).ToList();
            foreach (var l in l0)
            {
                yield return l;
            }

            var gi = sff.Aggregate(f.One, (acc, a) => acc * a.g.Pow(a.i));
            var c = f / gi;
            var cf = DeflateP(c);
            if (cf.q != 0)
            {
                foreach (var l in GianniTrager(cf.c, cf.q * q))
                {
                    yield return l;
                }
            }
        }
    }

    static void CheckSeparability<K>(KPoly<K> f, (KPoly<K> g, int q, int m)[] fSep)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = f.X;
        var fSep2 = fSep.Select(e => (gxp: e.g.Substitute(x.Pow(e.q)), e.q, e.m)).ToArray();
        var f0 = fSep2.Aggregate(x.One, (acc, e) => acc * e.gxp.Pow(e.m));
        if (f.Equals(f0))
        {
            Console.WriteLine("Prop (S1) pass");
        }
        else
        {
            Console.WriteLine("Prop (S1) fail");
            return;
        }

        var tuples = fSep2.Grid2D(fSep2).Where(e => !e.t1.gxp.Equals(e.t2.gxp)).ToArray();
        if (tuples.All(e => Ring.Gcd(e.t1.gxp, e.t2.gxp).Monic.Equals(f.One)))
        {
            Console.WriteLine("Prop (S2) pass");
        }
        else
        {
            Console.WriteLine("Prop (S2) fail");
            var pb = tuples.First(e => !Ring.Gcd(e.t1.gxp, e.t2.gxp).Monic.Equals(f.One));
            Console.WriteLine(pb);
            return;
        }

        var p = f.P;
        if (p == 0 || fSep2.All(e => e.m % p != 0))
        {
            Console.WriteLine("Prop (S3) pass");
        }
        else
        {
            Console.WriteLine("Prop (S3) fail");
            return;
        }

        if (fSep.All(e => e.g.Degree >= 1 && !Ring.Discriminant(e.g).IsZero()))
        {
            Console.WriteLine("Prop (S4) pass");
        }
        else
        {
            Console.WriteLine("Prop (S4) fail");
            var pb = fSep.First(e => e.g.Degree == 0 || Ring.Discriminant(e.g).IsZero());
            Console.WriteLine(pb);
            Console.WriteLine(Ring.Discriminant(pb.g));
            return;
        }

        if (tuples.All(e => !(e.t1.q, e.t1.m).Equals((e.t2.q, e.t2.m))))
        {
            Console.WriteLine("Prop (S5) pass");
        }
        else
        {
            Console.WriteLine("Prop (S5) fail");
            return;
        }

        Console.WriteLine("Successful Factorization");
        Console.WriteLine();
    }

    public static void SquareFreeFactorizationQ()
    {
        Monom.Display = MonomDisplay.Superscript;
        var x = FG.QPoly();

        {
            var f = (x + 1) * (x + 2) * (x + 3).Pow(2) * (x + 4).Pow(2) * (x + 5).Pow(3) * (x + 6).Pow(3);
            Console.WriteLine(f);
            Console.WriteLine("Musser Algo");
            var sff1 = MusserSFF(f).ToArray();
            Console.WriteLine(sff1.Glue("\n"));
            CheckSeparability(f, sff1);

            Console.WriteLine("Yun Algo");
            var sff2 = YunSFF(f).ToArray();
            Console.WriteLine(sff2.Glue("\n"));
            CheckSeparability(f, sff2);
        }

        {
            var f = x.Pow(3) * (x + 2).Pow(4) * (x.Pow(2) + 2 * x - 2).Pow(2) * (x.Pow(3) + 5);
            Console.WriteLine(f);
            Console.WriteLine("Musser Algo");
            var sff1 = MusserSFF(f).ToArray();
            Console.WriteLine(sff1.Glue("\n"));
            CheckSeparability(f, sff1);

            Console.WriteLine("Yun Algo");
            var sff2 = YunSFF(f).ToArray();
            Console.WriteLine(sff2.Glue("\n"));
            CheckSeparability(f, sff2);
        }
    }

    public static void SeparableFactorizationFp()
    {
        Monom.Display = MonomDisplay.Superscript;

        {
            var x = FG.ZPoly(3);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x + 2).Pow(4);
            Console.WriteLine(f);
            var l = GianniTrager(f, 1).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
        }

        {
            var x = FG.ZPoly(2);
            var f = x.Pow(2) * (x + 1).Pow(3) * (x.Pow(2) + 1).Pow(4);
            Console.WriteLine(f);
            var l = GianniTrager(f, 1).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
        }

        {
            var (x, t) = FG.FpT_Poly(3);
            // F = (X + 2T)7 (X3 + 2T)3 (X6 + T)
            var f = (x + 2 * t).Pow(7) * (x.Pow(3) + 2 * t).Pow(3) * (x.Pow(6) + t);
            Console.WriteLine(f);
            var l = GianniTrager(f, 1).ToArray();
            Console.WriteLine(l.Glue("\n"));
            CheckSeparability(f, l);
        }
    }
}