using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;

namespace FastGoat.UserGroup.Polynoms;

public static partial class IntFactorisation
{
    public static List<(KPoly<K> g, int q, int i)> MusserSFF<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var L = new List<(KPoly<K>, int, int)>();
        var c = Ring.FastGCD(f, f.Derivative);
        var i = 1;
        var g = f / c;
        while (g.Degree >= 1)
        {
            var p = Ring.FastGCD(c, g);
            c = c / p;
            if (g.Degree > p.Degree)
                L.Add(((g / p).Monic, 1, i));

            g = p;
            ++i;
        }

        return L;
    }
    
    public static List<(KPoly<K> g, int q, int i)> YunSFF<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var L = new List<(KPoly<K>, int, int)>() { (f.LT * f.One, 1, 1) };
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

    static (KPoly<K> c, int p) DeflateP<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var p = f.KZero.P;
        if (p == 0)
            throw new ArgumentException();

        var d = f.Degree;
        var decomp = IntExt.PrimesDecomposition(d).ToArray();
        if (!decomp.Contains(p))
            return (f, 1);

        var q = decomp.Count(i => i == p);
        var dq = (d + 1).Range();
        if (dq.Where(j => j % p != 0).All(j => f[j].IsZero()))
        {
            var coefs = dq.Where(j => j % p == 0).Select(j => f[j]).ToArray();
            return (new KPoly<K>(f.x, f.KZero, coefs), p);
        }

        return (f, 1);
    }

    // Page 334, Book ACFE, 18.4 Factorisation séparable
    public static IEnumerable<(KPoly<K> g, int q, int m)> GianniTrager<K>(KPoly<K> f, int q = 1)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (f.Degree != 0)
        {
            var sff = IntFactorisation.MusserSFF(f);
            var l0 = sff.Select(e => (e.g, q, m: e.i)).ToList();
            foreach (var l in l0)
            {
                yield return l;
            }

            var gi = sff.Aggregate(f.One, (acc, a) => acc * a.g.Pow(a.i));
            var c = f / gi;
            var cf = DeflateP(c);
            if (cf.p != 1)
            {
                foreach (var l in GianniTrager(cf.c, cf.p * q))
                {
                    yield return l;
                }
            }
        }
    }

    static EPoly<K>[] CanonicalBase<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = new EPoly<K>(f);
        return f.Degree.Range().Select(i => x.Pow(i)).ToArray();
    }

    static KMatrix<K> BerlekampMatrix<K>(KPoly<K> f, EPoly<K>[] baseCan, int q)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = baseCan.Length;
        var M = new K[n, n];
        var polys = baseCan.Select(g => g.Pow(q) - g).ToArray();
        foreach (var (i, j) in n.Range().Grid2D(n.Range()))
        {
            M[i, j] = polys[j][i];
        }

        return new(M);
    }

    public static KPoly<K>[] FrobeniusKernel<K>(KPoly<K> f, int q)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var baseCan = CanonicalBase(f);
        var bm = BerlekampMatrix(f, baseCan, q);
        var (nt, ns) = bm.NullSpace();
        var (m, n) = ns.Dim;
        var polys = n.Range().Select(j => m.Range().Select(i => ns[i, j] * baseCan[i]).Aggregate((a, b) => a + b).Poly)
            .ToArray();
        return polys;
    }

    static IEnumerable<KPoly<K>> FirrInternal<K>(KPoly<K> f, List<K> allF) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var polys = FrobeniusKernel(f, allF.Count);
        if (polys.Length > 1)
        {
            foreach (var (g, a) in polys.Where(g => g.Degree > 0).Grid2D(allF.Skip(1)))
            {
                var g_a = g - a;
                var gcd = Ring.FastGCD(f, g_a).Monic;
                if (gcd.Degree != 0)
                {
                    foreach (var f1 in FirrInternal(gcd, allF))
                        yield return f1;

                    foreach (var f1 in FirrInternal(f / gcd, allF))
                        yield return f1;

                    break;
                }
            }
        }
        else
        {
            yield return f;
        }
    }

    public static IEnumerable<KPoly<K>> Firr<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var acc = a0.One;
        var allF = new List<K>() { a0.Zero };
        do
        {
            allF.Add(acc);
            acc *= a0;
        } while (!acc.Equals(a0.One));

        return FirrInternal(f, allF);
    }

    public static List<(KPoly<K> g, int q, int m)> FirrFsep<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        List<(KPoly<K> g, int q, int m)> all = new();
        foreach (var (g, q, m) in GianniTrager(f))
            all.AddRange(Firr(g, a0).Select(g0 => (g0, q, m)));

        return all;
    }


    static EPoly<K> Mk<K>(EPoly<K> g, int k, int q) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (Int32.IsPow2(q))
        {
            var d = Int32.Log2(q);
            var g0 = g.Zero;
            var g1 = g;
            for (int i = 0; i < d * k; i++)
            {
                g0 += g1;
                g1 *= g1;
            }

            return g0;
        }
        else
        {
            var g0 = g.One;
            var g1 = g;
            for (int i = 0; i < k; i++)
            {
                g0 *= g1;
                g1 = g1.Pow(q);
            }

            return g0.Pow((q - 1) / 2) - 1;
        }
    }

    // A Computational Introduction to Number Theory and Algebra
    // Victor Shoup
    // 20.5 Factoring polynomials: Berlekamp’s algorithm page 541
    public static IEnumerable<KPoly<K>> BerlekampProbabilisticVShoup<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var acc = a0.One;
        var allF = new List<K>() { a0.Zero };
        var q = 0;
        do
        {
            allF.Add(acc);
            acc *= a0;
            ++q;
        } while (!acc.Equals(a0.One));

        var polys = FrobeniusKernel(f, q + 1);
        var r = polys.Length;

        var H0 = new List<KPoly<K>>() { f };
        var H1 = new List<KPoly<K>>();
        while (H0.Count < r)
        {
            var g = polys.Aggregate(f.Zero, (sum, gi) => sum + allF[IntExt.Rng.Next(q + 1)] * gi);
            H1.Clear();
            foreach (var h in H0)
            {
                var b = new EPoly<K>(h, g);
                var d = Ring.Gcd(Mk(b, 1, q + 1).Poly, h).Monic;
                if (d.Degree == 0 || d.Degree == h.Degree)
                    H1.Add(h.Monic);
                else
                    H1.AddRange(new[] { d, (h / d).Monic });
            }

            H0.Clear();
            H0.AddRange(H1);
        }

        return H0;
    }

    // AECF Algorithme de Berlekamp 353
    public static IEnumerable<KPoly<K>> BerlekampProbabilisticAECF<K>(KPoly<K> f, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        // var irrs = PolynomialFactorization.Firr(f, a0).Order().ToArray();
        // var r = irrs.Length.Range();
        // foreach (var (i, j) in r.Grid2D(r))
        // {
        //     var w = FG.EPoly(irrs[j], 'w');
        //     var df = f.Derivative.Substitute(w);
        //     var b0 = f / irrs[i] * irrs[i].Derivative;
        //     var bi = b0.Substitute(w) / df;
        //     Console.WriteLine(((i, j), bi)); // (Bi mod Fj) = (i==j ? 1 : 0)
        // }

        var acc = a0.One;
        var allF = new List<K>() { a0.Zero };
        do
        {
            allF.Add(acc);
            acc *= a0;
        } while (!acc.Equals(a0.One));

        var Gi = FrobeniusKernel(f, allF.Count);
        return BerlekampRec(f, Gi, allF);
    }

    // AECF Algorithme de Berlekamp 351
    static IEnumerable<KPoly<K>> BerlekampRec<K>(KPoly<K> F, KPoly<K>[] Gi, List<K> allF)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (Gi.Length == 1)
        {
            yield return F.Monic;
            yield break;
        }

        KPoly<K> H1;
        var q = allF.Count;
        var degF = F.Degree;
        while (true)
        {
            var G = Gi.Aggregate(F.Zero, (sum, gi) => sum + allF[IntExt.Rng.Next(q)] * gi);
            H1 = Ring.Gcd(F, G);
            if (H1.Degree > 0 && H1.Degree != degF)
                break;

            var a = new EPoly<K>(F, G);
            var H = Mk(a, 1, q).Poly;
            H1 = Ring.Gcd(F, H);
            if (H1.Degree > 0 && H1.Degree != degF)
                break;
        }

        var H2 = F / H1;

        var Gi1 = new HashSet<KPoly<K>>() { F.One };
        var Gi2 = new HashSet<KPoly<K>>() { F.One };

        var df = F.Derivative;
        var dh1h2 = new EPoly<K>(H1, H1.Derivative * H2).Inv();
        var h1dh2 = new EPoly<K>(H2, H1 * H2.Derivative).Inv();
        foreach (var gi in Gi)
        {
            var gidf = gi * df;

            var gi1 = (gidf * dh1h2).Poly;
            if (!gi1.IsZero())
                Gi1.Add(gi1.Monic);

            var gi2 = (gidf * h1dh2).Poly;
            if (!gi2.IsZero())
                Gi2.Add(gi2.Monic);
        }

        foreach (var f in BerlekampRec(H1, Gi1.ToArray(), allF))
            yield return f;

        foreach (var f in BerlekampRec(H2, Gi2.ToArray(), allF))
            yield return f;
    }

    static (int i, KPoly<K> Ei)[] DDF<K>(KPoly<K> F, int q) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var n = F.Degree;
        var X = FG.EPoly(F);
        var LX0 = new List<EPoly<K>>() { X };
        var f = F;
        for (int i = 1; i <= n; i++)
        {
            var Xi0 = LX0.Last();
            LX0.Add(Xi0.Pow(q));
        }

        var L = new List<(int i, KPoly<K> Ei)>();
        foreach (var (poly, i) in LX0.Select((ei, i) => (ei.Poly, i)).Skip(1))
        {
            var Ei = Ring.Gcd(poly - F.X, f).Monic;
            if (!Ei.Equals(f.One))
                L.Add((i, Ei));

            f = f / Ei;
        }

        return L.ToArray();
    }

    static IEnumerable<KPoly<K>> EDF<K>(KPoly<K> F, List<K> allF, int i) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (i == F.Degree)
        {
            yield return F;
            yield break;
        }

        var q = allF.Count;

        KPoly<K> H1;
        while (true)
        {
            var G = new KPoly<K>(F.x, F.KZero, F.Degree.Range().Select(i0 => allF[IntExt.Rng.Next(q)]).TrimSeq().ToArray());
            if (G.Degree == 0)
                continue;

            H1 = Ring.Gcd(F, G).Monic;
            if (H1.Degree > 0 && H1.Degree != F.Degree)
                break;

            var x = FG.EPoly(F);
            var H = Mk(G.Substitute(x), i, q).Poly.Monic;
            if (H.Degree == 0)
                continue;

            H1 = Ring.Gcd(F, H).Monic;
            if (H1.Degree > 0 && H1.Degree != F.Degree)
                break;
        }

        foreach (var fi in EDF(H1, allF, i))
            yield return fi;

        var H2 = (F / H1).Monic;
        foreach (var fi in EDF(H2, allF, i))
            yield return fi;
    }

    public static IEnumerable<KPoly<K>> CantorZassenhausVShoup<K>(KPoly<K> F, K a0)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var acc = a0.One;
        var allF = new List<K>() { a0.Zero };
        do
        {
            allF.Add(acc);
            acc *= a0;
        } while (!acc.Equals(a0.One));

        var q = allF.Count;
        foreach (var (i, ei) in DDF(F, q))
        {
            var r = ei.Degree / i;
            var H = new List<KPoly<K>>() { ei };
            while (H.Count < r)
            {
                var Hp = new HashSet<KPoly<K>>();
                foreach (var h in H)
                {
                    var a = new EPoly<K>(h,
                        new KPoly<K>(h.x, h.KZero, h.Degree.Range().Select(i0 => allF[IntExt.Rng.Next(q)]).TrimSeq().ToArray()));
                    var d = Ring.Gcd(Mk(a, i, q).Poly, h).Monic;
                    if (d.Equals(d.One) || d.Equals(h))
                        Hp.Add(h);
                    else
                        Hp.UnionWith(new[] { d, (h / d).Monic });
                }

                H = Hp.ToList();
            }

            foreach (var h in H)
                yield return h;
        }
    }

    public static IEnumerable<KPoly<K>> CantorZassenhausAECF<K>(KPoly<K> F, K a0) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var acc = a0.One;
        var allF = new List<K>() { a0.Zero };
        do
        {
            allF.Add(acc);
            acc *= a0;
        } while (!acc.Equals(a0.One));

        var q = allF.Count;
        foreach (var (i, ei) in DDF(F, q))
        {
            foreach (var fi in EDF(ei, allF, i))
                yield return fi;
        }
    }
}