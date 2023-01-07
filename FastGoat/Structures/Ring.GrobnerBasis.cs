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
            var Spol = Ring.SPolynomial(fi, fj);
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
        var ti = a.Indeterminates.Last();
        var da = a.DegreeOf(ti);
        var db = b.DegreeOf(ti);
        if (da != 0 && db != 0)
            throw new ArgumentException();

        var mnm = new Monom<Xi>(a.Indeterminates, ti, 1);
        var t = new Polynomial<K, Xi>(mnm, a.KOne);
        var gb = ReducedGrobnerBasis(t * a, (1 - t) * b);
        return gb.Last().Monic();
    }

}