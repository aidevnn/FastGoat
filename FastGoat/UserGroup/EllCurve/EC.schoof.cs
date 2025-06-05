using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;
using GFelt = FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>;

namespace FastGoat.UserGroup.EllCurve;

public static partial class EC
{
    public static (EllPoly<K> X, EllPoly<K> Y, EllPoly<K> Z) EllPoly<K>(K scalar, MonomOrder order = MonomOrder.Lex)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var s = new EllPoly<K>(scalar, order);
        return (s.X3, s.X2, s.X1);
    }

    public static EllPoly<Rational> Primitive(this EllPoly<Rational> f)
    {
        if (f.IsZero())
            return f;

        var arrGcd = f.Coefs.Values.Where(e => !e.IsZero()).Select(e => e.Absolute.Num).Distinct().Order().ToArray();
        var arrLcm = f.Coefs.Values.Select(e => e.Absolute.Denom).Distinct().Order().ToArray();
        return f * new Rational(f.LeadingDetails.lc.Sign * IntExt.LcmBigInt(arrLcm), IntExt.GcdBigInt(arrGcd));
    }

    public static EllPoly<ZnInt> ToZnInt(this EllPoly<Rational> f, int p)
    {
        var z = new ZnInt(p, 0);
        var coefs = f.Coefs.Select(e => (e.Key, e.Value.ToZnInt(p)))
            .Where(e => !e.Item2.IsZero()).ToDictionary(e => e.Key, e => e.Item2);
        return new(f.IndTriVar, z, coefs);
    }

    public static EllPoly<ZnBigInt> ToZnBigInt(this EllPoly<Rational> f, BigInteger p)
    {
        var z = new ZnBigInt(p, 0);
        var coefs = f.Coefs.Select(e => (e.Key, e.Value.ToZnBigInt(p)))
            .Where(e => !e.Item2.IsZero()).ToDictionary(e => e.Key, e => e.Item2);
        return new(f.IndTriVar, z, coefs);
    }

    public static (EllFracPoly<K> Y, EllFracPoly<K> X) EllFracPolyYX<K>((K a1, K a2, K a3, K a4, K a6) coefs,
        EllPoly<K> divPol)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var ind = divPol.IndTriVar;
        var (a1, a2, a3, a4, a6) = coefs;
        var x = new EllPoly<K>(ind, a1.Zero).X1;
        var y = x.X2;
        var lhs = y * y + a1 * x * y + a3 * y;
        var rhs = x.Pow(3) + a2 * x * x + a4 * x + a6;
        var (eq, sd) = (lhs - rhs, 2 * y + a1 * x + a3);
        var (Y, X) = new[] { eq.X2, eq.X1 }
            .Select(xi => new EllFracPoly<K>((eq, sd, divPol), xi, xi.One))
            .Deconstruct();
        return (Y, X);
    }

    public static (EllFracPoly<K> Y, EllFracPoly<K> X) EllFracPolyYX<K>(this EllGroup<K> ell, EllPoly<K> divPol)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EllFracPolyYX(ell.Coefs, divPol);
    }

    public static (EllFracPoly<K> Y, EllFracPoly<K> X) EllFracPolyYX<K>(this EllGroup<K> ell)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return ell.EllFracPolyYX(EllPoly(ell.a1).X.Zero);
    }

    public static (EllFracPoly<K> Y, EllFracPoly<K> X) EllFracPolyYX<K>(this EllCoefs<K> ell, EllPoly<K> divPol)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EllFracPolyYX(ell.Model, divPol);
    }

    public static (EllFracPoly<K> Y, EllFracPoly<K> X) EllFracPolyYX<K>(this EllCoefs<K> ell)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return ell.EllFracPolyYX(EllPoly(ell.a1).X.Zero);
    }

    public static (EllFracPoly<K> X, EllFracPoly<K> Y) EllFracPolyYX<K>(EllPoly<K> eqEll, EllPoly<K> sd, EllPoly<K> dvp)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (Y, X) = new[] { eqEll.X2, eqEll.X1 }
            .Select(xi => new EllFracPoly<K>((eqEll, sd, dvp), xi, xi.One))
            .Deconstruct();
        return (Y, X);
    }

    public static EllGroupSymb<K> ToEllGroupSymb<K>(this EllCoefs<K> E)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = new EllPoly<K>(E.a1).X1;
        return new(E, x.Zero);
    }

    public static long BSGSlong<T>(IGroup<T> g, T a, T b, double ord) where T : struct, IElt<T>
    {
        var (m, tmp1) = ((long)Double.Sqrt(ord) + 1, a); // TODO: faster BSGS
        var L = new Dictionary<T, long>() { [a] = 1 };
        for (long i = 1; i < m; i++)
        {
            if (tmp1.Equals(b))
                return i;

            tmp1 = g.Op(tmp1, a);
            L[tmp1] = i + 1;
        }

        if (tmp1.Equals(b))
            return m;

        var (c, tmp2) = (g.Invert(tmp1), b);
        for (long j = 1; j < m; j++)
        {
            tmp2 = g.Op(tmp2, c);
            if (L.TryGetValue(tmp2, out long i))
                return j * m + i;
        }

        throw new($"{g.Name} ord={ord}; a={a} b={b}");
    }

    public static int EllApBSGS(EllGroup<Rational> E, BigInteger p, LogLevel log = LogLevel.Off)
    {
        if (p < 200 || E.Disc.ToZnBigInt(p).IsZero()) // TODO: #Ep when disc = 0 mod p
            return EllAp(E, (int)p);

        var d = double.Sqrt((double)p);
        var (pmin, pmax) = ((double)p + 1 - 2 * d, (double)p + 1 + 2 * d);

        (BigInteger L, BigInteger N) = (1, -1);
        (BigInteger L_, BigInteger N_) = (1, -1);
        var Ep = E.ToZnBigInt(p);
        var Es = Ep.ToShortWeierstrassForm();

        var g = 1000.SeqLazy().Select(_ => DistributionExt.Dice(BigInteger.One * 2, p - 1))
            .Where(g => LegendreJacobiBigint(g, p) != 1)
            .Select(g => new ZnBigInt(p, g))
            .Where(g => !(4 * (Es.a4 * g.Pow(2)).Pow(3) + 27 * (Es.a6 * g.Pow(3)).Pow(2)).IsZero())
            .First();

        // Mestre's Theorem, quadratic twist of Ep
        var Ep_ = new EllGroup<ZnBigInt>(Es.a4 * g.Pow(2), Es.a6 * g.Pow(3));
        int ct = 0, nbCands = 0, nbCands_ = 0;
        var set = new HashSet<BigInteger>();
        var set_ = new HashSet<BigInteger>();
        var cands = new List<BigInteger>();
        var cands_ = new List<BigInteger>();
        while (ct < 2 * d)
        {
            ++ct;
            var P = RandPt(Ep, set);
            var ordP = BSGSlong(Ep, P, Ep.O, pmax);

            L = LcmBigInt(ordP, L);
            cands = CandidatsN(pmin, pmax, L);
            if (cands.Count >= 1)
                N = cands[0] / L;

            var NL = N * L;
            var arrCands = $"[{cands.Glue(", ")}]";
            nbCands = cands.Count;

            var P_ = RandPt(Ep_, set_);
            var ordP_ = BSGSlong(Ep_, P_, Ep_.O, pmax);
            L_ = LcmBigInt(ordP_, L_);
            cands_ = CandidatsN(pmin, pmax, L_);
            if (cands_.Count >= 1)
                N_ = cands_[0] / L_;

            var NL_ = N_ * L_;
            var arrCands_ = $"[{cands_.Glue(", ")}]";
            nbCands_ = cands_.Count;

            if (log == LogLevel.Level2)
            {
                Console.WriteLine(new { p, d, pmin, pmax });
                Console.WriteLine(new { NL, N, L, arrCands, nbCands });
                Console.WriteLine(new { NL_, N_, L_, arrCands_, nbCands_ });
            }

            if (N == 1 || N_ == 1)
                break;

            if (nbCands == 1 && L % N == 0 && (p - 1) % N == 0 && cands_.Contains(2 * p + 2 - NL))
            {
                var N0 = (2 * p + 2 - NL) / L_;
                if (L_ % N0 == 0 && (p - 1) % N0 == 0)
                    break;
            }

            if (nbCands_ == 1 && L_ % N_ == 0 && (p - 1) % N_ == 0 && cands.Contains(2 * p + 2 - NL_))
            {
                var N0_ = (2 * p + 2 - NL_) / L;
                if (L % N0_ == 0 && (p - 1) % N0_ == 0)
                    break;
            }
        }

        int caseNum, outputNum;
        if (N == 1)
            (caseNum, outputNum) = (1, 1);
        else if (N_ == 1)
            (caseNum, outputNum) = (2, 2);
        else if (L % N == 0 && (p - 1) % N == 0 && nbCands == 1 && L_ < L)
            (caseNum, outputNum) = (3, 1);
        else if (L_ % N_ == 0 && (p - 1) % N_ == 0 && nbCands_ == 1 && L < L_)
            (caseNum, outputNum) = (4, 2);
        else if (nbCands == 1 && nbCands_ == 1 && 2 * (p + 1) == cands[0] + cands_[0])
            (caseNum, outputNum) = (5, 3);
        else if (nbCands == 2 && nbCands_ == 2)
        {
            if (L_ < L)
                (caseNum, outputNum) = (6, 1);
            else if (L < L_)
                (caseNum, outputNum) = (7, 2);
            else
                (caseNum, outputNum) = (8, 2);
        }
        else
            (caseNum, outputNum) = (9, 4);

        if (outputNum == 1 || outputNum == 3)
        {
            var Ap = (int)cands.Where(NL => (p - 1) % (NL / L) == 0).Select(NL => p + 1 - NL).First();
            if (log != LogLevel.Off)
                Console.WriteLine($"Ep = {Ep} Ap = {Ap} Case = [{caseNum}]#{outputNum}");

            return Ap;
        }
        else if (outputNum == 2)
        {
            var Ap = (int)cands_.Where(NL => (p - 1) % (NL / L_) == 0).Select(NL => -(p + 1 - NL)).First();
            if (log != LogLevel.Off)
                Console.WriteLine($"Ep = {Ep} Ap = {Ap} Case = [{caseNum}]#2");

            return Ap;
        }
        else
            throw new();
    }

    public static (Dictionary<int, EllFracPoly<K>> psi, Dictionary<int, EllFracPoly<K>> f)
        DivisionPolynomial<K>(EllGroup<K> E, int nmax) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (Y, X) = E.EllFracPolyYX();
        var (X2, X3, X4, X5, X6) = 5.SeqLazy(2).Select(i => X.Pow(i)).Deconstruct();
        var eCoefs = E.ToEllCoefs().ToEllCoefs(X);
        var (a1, a2, a3, a4, a6) = eCoefs.Model;
        var (b2, b4, b6, b8) = eCoefs.B_Invariants;

        var psi = new Dictionary<int, EllFracPoly<K>>();
        (psi[0], psi[1], psi[2]) = (Y.Zero, Y.One, 2 * Y + a1 * X + a3);
        psi[3] = 3 * X4 + b2 * X3 + 3 * b4 * X2 + 3 * b6 * X + b8;
        psi[4] = psi[2] * (2 * X6 + b2 * X5 + 5 * b4 * X4 + 10 * b6 * X3 + 10 * b8 * X2 + (b2 * b8 - b4 * b6) * X +
                           (b4 * b8 - b6 * b6));

        for (int n = 2; n <= nmax / 2; n++)
        {
            psi[2 * n] = psi[n] * (psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2)) / psi[2];
            if (2 * n + 1 <= nmax)
                psi[2 * n + 1] = psi[n + 2] * psi[n].Pow(3) - psi[n + 1].Pow(3) * psi[n - 1];
        }

        var f = psi.ToDictionary(e => e.Key, e => e.Key % 2 == 0 ? (e.Value / psi[2]) : e.Value);
        return (psi, f);
    }

    public static EllPt<EllFracPoly<K>> NPt<K>(int n, Dictionary<int, EllFracPoly<K>> psi)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var psi2 = psi[2];
        var X = EllFracPolyYX(psi2.EqEll, psi2.SD, psi2.DivPol).X;
        var (a1, a3) = (psi2.Num[(0, 0, 1)], psi2.Num.ConstTerm);

        var nPX = X - psi[n - 1] * psi[n + 1] / psi[n].Pow(2);

        var num = psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2);
        var denom = 2 * psi[2] * psi[n].Pow(3);
        var nPY = num / denom - (a1 * nPX + a3) / 2;

        return new(nPX, nPY);
    }

    public static GFelt[] Roots(KPoly<GFelt> P)
    {
        var a = P.KOne.X;
        var q = BigInteger.Pow(P.P, a.F.Degree);
        var facts = IntFactorisation.MusserSFF(P)
            .SelectMany(f => IntFactorisation.CantorZassenhausAECF(f.g, a, q))
            .ToArray();

        if (facts.All(f => f.Degree == 1))
            return facts.Select(f => -f[0]).ToArray();

        return [];
    }

    public static GFelt[] Roots(KPoly<GFelt> P, GFelt a)
    {
        var q = BigInteger.Pow(P.P, a.F.Degree);
        var facts = IntFactorisation.CantorZassenhausAECF(P, a, q).ToArray();

        if (facts.All(f => f.Degree == 1))
            return facts.Select(f => -f[0]).ToArray();

        return [];
    }

    static (GFelt[], GFelt g) NTorsSplittingFp(KPoly<ZnInt> P)
    {
        var p = P.P;
        var a0 = NumberTheory.PrimitiveRootMod(p) * P.KOne;
        var facts = IntFactorisation.MusserSFF(P)
            .SelectMany(f => IntFactorisation.CantorZassenhausAECF(f.g, a0, p))
            .ToArray();

        if (facts.Length == 0)
            throw new();

        var fact = facts.MaxBy(f => f.Degree);
        var (X, a) = FG.EPolyXc(fact, 'a');
        var g = NumberTheory.PrimitiveRoot(a);
        // Irreductibles factors of division polynomial in Fp
        // when f = f1 * f2 * ... * fi orderer by degree
        // it seems that // deg(f1) | deg(f2) | ... | deg(fi) | deg(divPol)
        // and fi seems to be the minimal polynomial of the primitive element
        // of the splitting field of f in Fp[a]
        // TODO: proof
        var roots = facts.SelectMany(f => Roots(f.Substitute(X), g)).Distinct().ToArray();
        if (roots.Length == 0)
            throw new();

        return (roots, g);
    }

    static (HashSet<EllPt<GFelt>>, GFelt g)
        NTorsExtensionFp(EllGroup<Rational> Efq, KPoly<Rational> P, int p, LogLevel log = LogLevel.Level1)
    {
        var (roots0, g0) = NTorsSplittingFp(P.ToZnPoly(p));
        var Efq0 = Efq.ToGF(g0.X);
        var (a01, a02, a03, a04, a06) = Efq0.Coefs;
        var seq0 = roots0.Select(x => (x, b: a01 * x + a03, c: -(x.Pow(3) + a02 * x * x + a04 * x + a06)))
            .Select(e => (e.x, e.b, e.c, delta: e.b.Pow(2) - 4 * e.c))
            .Select(e => (e.x, e.b, e.c, e.delta, resQuad: NumberTheory.LegendreJacobiGf(e.delta))).ToList();
        var test0 = seq0.All(e => e.delta.IsZero() || e.resQuad.IsOne());
        if (test0)
        {
            if (log != LogLevel.Off)
            {
                Console.WriteLine($"x(P) and y(P) in GF({p}^{g0.F.Degree})");
                Console.WriteLine($"    MinPol   of GF({p}^{g0.F.Degree}) = {g0.F.SubstituteChar('X')}");
                Console.WriteLine($"    PrimRoot of GF({p}^{g0.F.Degree}) = {g0}");
            }

            var sqrts0 = seq0.Select(e => (e.x, e.b, e.c, e.delta, sqrtDelta: NumberTheory.SqrtFqANTV1(e.delta, g0)))
                .ToArray();
            return (sqrts0.SelectMany(e => new[]
            {
                new EllPt<GFelt>(e.x, (-e.b + e.sqrtDelta) / 2),
                new EllPt<GFelt>(e.x, (-e.b - e.sqrtDelta) / 2)
            }).ToHashSet(), g0);
        }
        else
        {
            if (log != LogLevel.Off)
            {
                Console.WriteLine($"x(P) in GF({p}^{g0.F.Degree})");
                Console.WriteLine($"    MinPol   of GF({p}^{g0.F.Degree}) = {g0.F.SubstituteChar('X')}");
                Console.WriteLine($"    PrimRoot of GF({p}^{g0.F.Degree}) = {g0}");
            }

            var delta = seq0.First(e => !e.delta.IsZero() && !e.resQuad.IsOne()).delta;
            var Y = FG.KPoly('Y', g0.X);
            var (r, a, _) = IntFactorisation.PrimitiveElt(Y.Pow(2) - delta);
            var a0 = FG.EPoly(r, delta.Poly.x);
            var a1 = a.Substitute(a0);
            var g1 = NumberTheory.PrimitiveRoot(a0.X);
            var seq1 = seq0.Select(e => (x: e.x.Substitute(a1), b: e.b.Substitute(a1), c: e.c.Substitute(a1),
                    delta: e.delta.Substitute(a1)))
                .Select(e => (e.x, e.b, e.c, e.delta, resQuad: NumberTheory.LegendreJacobiGf(e.delta))).ToList();

            var test1 = seq1.All(e => e.delta.IsZero() || e.resQuad.IsOne());
            if (!test1)
                throw new("#1");

            if (log != LogLevel.Off)
            {
                Console.WriteLine($"y(P) in GF({p}^{g1.F.Degree})");
                Console.WriteLine($"    MinPol   of GF({p}^{g1.F.Degree}) = {g1.F.SubstituteChar('X')}");
                Console.WriteLine($"    PrimRoot of GF({p}^{g1.F.Degree}) = {g1}");
            }

            var sqrts1 = seq1.Select(e => (e.x, e.b, e.c, e.delta, sqrtDelta: NumberTheory.SqrtFqANTV1(e.delta, g1)));
            return (sqrts1.SelectMany(e => new[]
                {
                    new EllPt<GFelt>(e.x, (-e.b + e.sqrtDelta) / 2),
                    new EllPt<GFelt>(e.x, (-e.b - e.sqrtDelta) / 2)
                })
                .ToHashSet(), g1);
        }
    }

    public static (ConcreteGroup<EllPt<GFelt>> nTors, GFelt g)
        EllFpNTors(EllGroup<Rational> Efq, KPoly<Rational> P, int p, int l, LogLevel log = LogLevel.Level1)
    {
        if (log != LogLevel.Off)
            GlobalStopWatch.AddLap();

        var (pts, g) = NTorsExtensionFp(Efq, P, p, log);
        var Efq1 = Efq.ToGF(g);
        if (pts.Any(pt => !Efq1.Contains(pt.X, pt.Y)))
            throw new("#2");

        var orders = pts.ToDictionary(pt => pt, pt => Group.Cycle(Efq1, pt));
        var e1 = orders.OrderByDescending(e => e.Key.X.Degree + e.Key.Y.Degree)
            .First(e => e.Value.Count == l).Key;
        var e2 = orders.OrderByDescending(e => e.Key.X.Degree + e.Key.Y.Degree)
            .First(e => e.Value.Count == l && !e.Value.Keys.Contains(e1)).Key;
        var nTors = Group.Generate(Efq1, [e1, e2]);
        nTors.Name = $"{l}-Tors({nTors.Name})";

        if (log != LogLevel.Off)
        {
            DisplayGroup.HeadGenerators(nTors);
            var abType = Group.AbelianGroupType(nTors);
            Console.WriteLine($"{nTors} ~ {abType.Glue(" x ", "C{0}")}");
            GlobalStopWatch.Show();
            Console.WriteLine();
        }

        return (nTors, g);
    }

    public static Dictionary<int, List<int>> SmallPrimesList(Rational N)
    {
        var nmax = 2 * double.Sqrt(N);
        var listL = new Dictionary<int, List<int>>();
        foreach (var p in Primes10000.Where(p => p <= nmax && !N.Mod(p).IsZero()))
        {
            var lmax = 4 * double.Sqrt(p);
            listL[p] = new();
            var L = 1;
            foreach (var l in Primes10000.Where(l => p % l != 0))
            {
                listL[p].Add(l);
                L *= l;
                if (L > lmax)
                    break;
            }
        }

        return listL;
    }

    public static EllPt<GFelt> Frob(EllPt<GFelt> pt, int n = 1)
    {
        if (pt.IsO)
            return pt;

        var p = pt.X.P;
        return new(pt.X.FastPow(BigInteger.Pow(p, n)), pt.Y.FastPow(BigInteger.Pow(p, n)));
    }

    public static int FrobTrace(EllGroup<GFelt> E, EllPt<GFelt>[] nTors, int l) =>
        l.SeqLazy().First(t => nTors.All(pt =>
        {
            if (pt.IsO)
                return true;

            var p = pt.X.P;
            var phi2 = Frob(Frob(pt));
            var t_phi = Frob(E.Times(pt, -t));
            var q = E.Times(pt, p);
            return E.Op(phi2, E.Op(t_phi, q)).IsO;
        }));

    public static (Dictionary<int, KPoly<Rational>> divPolys,
        Dictionary<int, Dictionary<int, (GFelt g, EllPt<GFelt> e1, EllPt<GFelt> e2)>> allBasis)
        EllApFrobTrace(BigInteger[] curve)
    {
        GlobalStopWatch.AddLap();
        var E = EllCoefs(curve);
        var N = EllTateAlgorithm(EllCoefs(curve)).N;
        var Ell = E.ToEllGroup();
        var allList = SmallPrimesList(N);
        var pmax = allList.Keys.Max();
        var lmax = allList.Max(e => e.Value.Max());
        Console.WriteLine($"{Ell}");

        var (psi0, fdiv0) = DivisionPolynomial(Ell, lmax + 3);
        var Pt2 = NPt(2, psi0);
        fdiv0[2] = new(Pt2.X.Reduction, Pt2.X.Denom, Pt2.X.Num.One);
        var divPolys = fdiv0.ToDictionary(e => e.Key, e => e.Value.Num.ToKPolyX1())
            .ToDictionary(e => e.Key, e => e.Value.PrimitiveZPoly());
        divPolys.Println("divPolys");

        Console.WriteLine($"N = {N} pmax = {pmax} listMax = {lmax}");
        var allBasis = new Dictionary<int, Dictionary<int, (GFelt g, EllPt<GFelt> e1, EllPt<GFelt> e2)>>();
        foreach (var (p, listL) in allList)
        {
            if (p < 5)
                continue;

            var ap = EllAp(Ell, p);
            var frobTr = new Dictionary<int, int>();
            var basis = new Dictionary<int, (GFelt g, EllPt<GFelt> e1, EllPt<GFelt> e2)>();
            foreach (var l in listL)
            {
                if (p % l == 0)
                    continue;

                var Ep = Ell.ToGF(p);
                Ep.Field = $"GF({p})";
                Console.WriteLine($"{l}-torsion of {Ep} from {Ell}");
                Console.WriteLine($"{l}-divPol mod {p} = {divPolys[l].ToZnPoly(p)}");
                var (nTors, g) = EllFpNTors(Ell, divPolys[l], p, l);
                var Egf = Ell.ToGF(g);
                var (e1, e2) = nTors.GetGenerators().Deconstruct();
                frobTr[l] = FrobTrace(Egf, nTors.GetGenerators().ToArray(), l);
                basis[l] = (g, e1, e2);
            }

            allBasis[p] = basis;
            frobTr.Println("Frob Traces");
            var keys = frobTr.Keys.ToArray();
            var values = keys.Select(k => frobTr[k]).ToArray();
            var crtTable = NumberTheory.CrtTable(keys);
            var L = keys.Aggregate((li, lj) => li * lj);
            var crt = NumberTheory.CRT(values, crtTable, L);
            var ap1 = crt < L / 2 ? crt : crt - L;
            Console.WriteLine($"p = {p} ap = {ap} crt = {crt} ap1 = {ap1} L = {L} Check:{ap == ap1}");
            if (ap != ap1)
                throw new();

            Console.WriteLine();
        }

        GlobalStopWatch.Show($"End FrobTrace {Ell}");
        Console.WriteLine();
        return (divPolys, allBasis);
    }

    public static void EllApSchoof(BigInteger[] curve)
    {
        var E = EllCoefs(curve);
        var El = E.ToLongWeierstrassForm();
        var Ell = El.ToEllGroup();

        GlobalStopWatch.AddLap();
        var N = EllTateAlgorithm(E).N;
        var allList = SmallPrimesList(N);
        var pmax = allList.Keys.Max();
        var lmax = allList.Max(e => e.Value.Max());
        Console.WriteLine($"{E.ToEllGroup()} => {Ell} Conductor N = {N} j-Inv = {E.J_Invariant}");

        var (psi0, fdiv0) = DivisionPolynomial(Ell, lmax + 3);
        var Pt2 = NPt(2, psi0);
        var x = Pt2.X.Num.X1;
        fdiv0[2] = new(Pt2.X.Reduction, Pt2.X.Denom, Pt2.X.Num.One);
        var divPolys = fdiv0.ToDictionary(e => e.Key, e => e.Value.Num.Primitive());
        divPolys.Println("divPolys");

        Console.WriteLine($"N = {N} pmax = {pmax} listMax = {lmax}");
        var frobTr = new Dictionary<int, Dictionary<int, int>>();
        foreach (var (p, listL) in allList)
        {
            if (p <= 3)
                continue;

            var ap = EllAp(Ell, p);
            var pFrobTr = new Dictionary<int, int>();
            var g = NumberTheory.PrimitiveRootFp(p);
            var Ep = Ell.ToEllCoefs().ToZnInt(p);
            GlobalStopWatch.AddLap();
            foreach (var l in listL)
            {
                GlobalStopWatch.AddLap();
                var psi = divPolys[l].ToZnInt(p);

                var facts = IntFactorisation.FirrFsepCantorZassenhausAECF(psi.ToKPolyX1(), g, p);
                facts.Println($"psi = {psi}");
                var fPsi = facts.OrderDescending().First().g.Substitute(x.ToZnInt(p));
                // Irreductibles factors of division polynomial in Fp
                // when f = f1 * f2 * ... * fi orderer by degree
                // it seems that // deg(f1) | deg(f2) | ... | deg(fi) | deg(divPol)
                // and fi seems to be the minimal polynomial of the primitive element
                // of the splitting field of f in Fp[a]
                // TODO: proof
                var Erl = new EllGroupSymb<ZnInt>(Ep, fPsi);
                Erl.CheckValidity = false;
                Console.WriteLine($"p={p} l={l} {Erl}");
                Console.WriteLine($"{Erl.Eq}");
                Console.WriteLine($"psi  = {psi}");
                Console.WriteLine($"fPsi = {fPsi}");

                var pt = Erl.Pt;
                var p_Pt = Erl.Times(pt, p % l);
                var phi2 = Erl.FrobRl(pt, 2);
                var add_phi2_p = Erl.Op(phi2, p_Pt);

                var phi = Erl.Invert(Erl.FrobRl(pt));
                var t_phi = Erl.O;
                foreach (var t in l.SeqLazy())
                {
                    var eqFrob = Erl.Op(add_phi2_p, t_phi);
                    if (eqFrob.IsO)
                    {
                        pFrobTr[l] = t;
                        Console.WriteLine($"phi(P)^2 - t*phi(P) + p*P = O");
                        Console.WriteLine($"    t = {t}");
                        break;
                    }

                    t_phi = Erl.Op(phi, t_phi);
                }

                if (!pFrobTr.ContainsKey(l))
                    throw new();

                GlobalStopWatch.Show();
                Console.WriteLine();
            }

            pFrobTr.Println("Frob Traces");
            var keys = pFrobTr.Keys.ToArray();
            var values = keys.Select(k => pFrobTr[k]).ToArray();
            var crtTable = NumberTheory.CrtTable(keys);
            var L = keys.Aggregate((li, lj) => li * lj);
            var crt = NumberTheory.CRT(values, crtTable, L);
            var ap1 = crt < L / 2 ? crt : crt - L;
            Console.WriteLine($"p = {p} ap = {ap} crt = {crt} ap1 = {ap1} L = {L} Check:{ap == ap1}");
            GlobalStopWatch.Show($"EllApSchoof({Ep.ToEllGroupSymb()})");
            if (ap != ap1)
                throw new();

            Console.WriteLine();
        }

        GlobalStopWatch.Show($"End EllApSchoof {Ell}");
        Console.WriteLine();
    }
}