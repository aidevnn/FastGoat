using FastGoat.Commons;
using FastGoat.Structures.VecSpace;

namespace FastGoat.Structures;

public static partial class Ring
{
    static Polynomial<K, Xi> Reduction<K>(Polynomial<K, Xi> p, Polynomial<K, Xi>[] bs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var f = p;
        var allLT = bs.Select(e => (e, e.LeadingDetails)).ToArray();
        var stop = false;
        while (!f.IsZero() && !stop)
        {
            var l = f.LeadingDetails;
            stop = true;
            foreach (var (fk, lk) in allLT)
            {
                var (b, m) = l.lm.Div(lk.lm);
                if (b)
                {
                    var c = l.lc / lk.lc;
                    f -= new Polynomial<K, Xi>(m, c) * fk;
                    stop = false;
                    break;
                }
            }
        }

        return f;
    }

    public static Polynomial<K, Xi> TotalReduction<K>(Polynomial<K, Xi> p, Polynomial<K, Xi>[] bs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var f0 = p.Zero;
        var f = p;
        while (!f.IsZero())
        {
            f = Reduction(f, bs);
            if (f.IsZero())
                return f0;

            var lt = f.LeadingDetails.lt;
            f0 += lt;
            f -= lt;
        }

        return f0;
    }

    static IEnumerable<Polynomial<K, Xi>> MinGBasis<K>(IEnumerable<Polynomial<K, Xi>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (!f.Any())
            return Enumerable.Empty<Polynomial<K, Xi>>();

        var p = f.Max();
        var g1 = MinGBasis(f.Except(new[] { p })).ToList();
        var (lc, lm, _) = p.LeadingDetails;
        p /= lc;
        if (g1.All(e => !lm.Div(e.LeadingDetails.lm).Item1))
            g1.Add(p);

        return g1;
    }

    static IEnumerable<Polynomial<K, Xi>> ReducedBasis<K>(IEnumerable<Polynomial<K, Xi>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var g0 = MinGBasis(f).ToArray();
        return g0.Select(p => TotalReduction(p, g0.Except(new[] { p }).ToArray())).Where(p => !p.IsZero());
    }

    // Syzygie Polynomial
    public static (Polynomial<K, Xi> pi, Polynomial<K, Xi> pj, Polynomial<K, Xi> Sij) SPolynomial<K>(Polynomial<K, Xi> fi,
        Polynomial<K, Xi> fj)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (lci, lmi, _) = fi.LeadingDetails;
        var (lcj, lmj, _) = fj.LeadingDetails;
        var (szi, szj) = Monom<Xi>.Reduce(lmi, lmj);
        var pi = new Polynomial<K, Xi>(szi, lcj);
        var pj = new Polynomial<K, Xi>(szj, lci);
        var Sij = pi * fi - pj * fj;
        return (pi, pj, Sij);
    }

    // Jean-Charles Faugere version of Buchberger algorithm for Grobner Basis
    // Jean-Charles Faugere version of Buchberger algorithm for Groebner Basis
    static Polynomial<K, Xi>[] BuchbergerAlgorithm<K>(params Polynomial<K, Xi>[] f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var m = f.Length;
        var g = f.ToList();
        var Pij = m.Range().Grid2D(m.Range()).Where(e => e.t1 < e.t2).Select(e => (g[e.t1], g[e.t2])).ToList();
        while (Pij.Count != 0)
        {
            var (fi, fj) = Pij[0];
            Pij.RemoveAt(0);
            var Spol = SPolynomial(fi, fj);
            var fk = Reduction(Spol.Sij, g.ToArray());
            if (!fk.IsZero())
            {
                ++m;
                Pij.AddRange(g.Select(fh => (fh, fk)));
                g.Add(fk);
            }
        }

        return g.ToArray();
    }

    public static Polynomial<K, Xi>[] GrobnerBasis<K>(params Polynomial<K, Xi>[] f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return BuchbergerAlgorithm(f);
    }

    public static Polynomial<K, Xi>[] GroebnerBasis<K>(params Polynomial<K, Xi>[] f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return BuchbergerAlgorithm(f);
    }

    public static Polynomial<K, Xi>[] ReducedGrobnerBasis<K>(params Polynomial<K, Xi>[] f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return ReducedBasis(GrobnerBasis(f)).Order().ToArray();
    }

    public static Polynomial<K, Xi>[] ReducedGroebnerBasis<K>(params Polynomial<K, Xi>[] f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return ReducedBasis(GrobnerBasis(f)).Order().ToArray();
    }

    public static Polynomial<K, Xi> LcmPolynomial<K>(Polynomial<K, Xi> a, Polynomial<K, Xi> b)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        a.Indeterminates.ExtendPrepend(new Xi("_t_"));
        var ti = a.Indeterminates.First();
        var da = a.DegreeOf(ti);
        var db = b.DegreeOf(ti);
        if (da != 0 && db != 0)
            throw new ArgumentException();

        var t = ti.ToPolynomial(a.One);
        var ord = a.Indeterminates.Order;
        a.Indeterminates.SetOrder(MonomOrder.Lex);

        var gb = ReducedGrobnerBasis(t * a, (1 - t) * b);
        var lcm = gb.First(e => e.Coefs.Keys.All(xi => xi[ti] == 0));
        a.Indeterminates.SetOrder(ord);
        a.Indeterminates.Remove(ti);
        return lcm.Monic();
    }

    public static Polynomial<K, Xi> GcdPolynomial<K>(Polynomial<K, Xi> a, Polynomial<K, Xi> b)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return ((a * b) / LcmPolynomial(a, b)).Monic();
    }

    public static Polynomial<K, Xi> LcmPolynomials<K>(Polynomial<K, Xi> a, params Polynomial<K, Xi>[] ps)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (ps.Length == 0)
            return a;

        return LcmPolynomial(a, LcmPolynomials(ps[0], ps.Skip(1).ToArray()));
    }

    public static Polynomial<K, Xi> GcdPolynomials<K>(Polynomial<K, Xi> a, params Polynomial<K, Xi>[] ps)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (ps.Length == 0)
            return a;

        return GcdPolynomial(a, GcdPolynomials(ps[0], ps.Skip(1).ToArray()));
    }
}