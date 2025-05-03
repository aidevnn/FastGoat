using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Integers;
using static FastGoat.Commons.IntExt;
using GFelt = FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>;

namespace FastGoat.UserGroup.EllCurve;

public static class EC
{
    [Flags]
    public enum TorsionMeth
    {
        None = 0,
        Fp = 1,
        NagellLutz = 2,
        Both = 3
    }

    public static bool IsIntegral(this EllPt<Rational> pt) => pt.IsO || (pt.X.IsInteger() && pt.Y.IsInteger());

    public static BigInteger Height(this EllPt<Rational> pt)
    {
        if (pt.IsO)
            return 1;

        return BigInteger.Max(pt.X.Absolute.Num, pt.X.Absolute.Denom);
    }

    public static T[] CurveArray<T>(params T[] ts) where T : struct, IElt<T>, IRingElt<T>
    {
        if (ts.Length == 2)
        {
            var z = ts[0].Zero;
            return [z, z, z, ts[0], ts[1]];
        }
        else if (ts.Length == 3)
        {
            var z = ts[0].Zero;
            return [z, ts[0], z, ts[1], ts[2]];
        }
        else if (ts.Length == 5)
            return ts;
        else
            throw new();
    }

    public static EllGroup<Rational> EllGroup(BigInteger[] curve)
    {
        var (a1, a2, a3, a4, a6) = CurveArray(curve.Select(e => new Rational(e)).ToArray()).Deconstruct();
        return new(a1, a2, a3, a4, a6);
    }

    public static EllGroup<Rational> EllGroup(int[] curve)
    {
        var (a1, a2, a3, a4, a6) = CurveArray(curve.Select(e => new Rational(e)).ToArray()).Deconstruct();
        return new(a1, a2, a3, a4, a6);
    }

    public static EllCoefs<Rational> EllCoefs(BigInteger[] curve)
    {
        var (a1, a2, a3, a4, a6) = CurveArray(curve.Select(e => new Rational(e)).ToArray()).Deconstruct();
        return new(a1, a2, a3, a4, a6);
    }

    public static EllCoefs<Rational> EllCoefs(int[] curve)
    {
        var (a1, a2, a3, a4, a6) = CurveArray(curve.Select(e => new Rational(e)).ToArray()).Deconstruct();
        return new(a1, a2, a3, a4, a6);
    }

    public static ConcreteGroup<EllPt<ZnBigInt>> EllFp(EllGroup<ZnBigInt> E)
    {
        var (A, B) = (E.ShortForm.A, E.ShortForm.B);
        var p = (int)A.Mod;
        var ell = p.Range().Select(k => new ZnBigInt(p, k))
            .Select(x => (x, y2: x.Pow(3) + A * x + B))
            .Select(e => (e.x, y: NumberTheory.SqrtModANTV1(e.y2.K, p) * e.x.One))
            .Select(e => E.ConvertFromShort(new(e.x, e.y)))
            .Where(e => E.Contains(e.X, e.Y))
            .Order()
            .ToArray();

        return Group.Generate(E, ell);
    }

    public static EllGroup<ZnInt> ToZnInt(this EllGroup<Rational> E, int p)
    {
        var (a1, a2, a3, a4, a5) = E.Coefs;
        return new(a1.ToZnInt(p), a2.ToZnInt(p), a3.ToZnInt(p), a4.ToZnInt(p), a5.ToZnInt(p));
    }

    public static EllGroup<GFelt> ToGF(this EllGroup<Rational> E, BigInteger q, char a = 'a')
    {
        var (a1, a2, a3, a4, a5) = E.Coefs;
        return new(a1.ToGF(q, a), a2.ToGF(q, a), a3.ToGF(q, a), a4.ToGF(q, a), a5.ToGF(q, a));
    }

    public static EllCoefs<GFelt> ToGF(this EllCoefs<Rational> E, BigInteger q, char a = 'a')
    {
        var (a1, a2, a3, a4, a5) = E.Model;
        return new(a1.ToGF(q, a), a2.ToGF(q, a), a3.ToGF(q, a), a4.ToGF(q, a), a5.ToGF(q, a));
    }

    public static EllGroup<GFelt> ToGF(this EllGroup<GFelt> E, BigInteger q, char a = 'a')
    {
        var (a1, a2, a3, a4, a6) = E.Coefs;
        return new(a1.ToGF(q, a), a2.ToGF(q, a), a3.ToGF(q, a), a4.ToGF(q, a), a6.ToGF(q, a));
    }

    public static EllGroup<GFelt> ToGF(this EllGroup<GFelt> E, GFelt a)
    {
        var (a1, a2, a3, a4, a6) = E.Coefs;
        return new(a1.Substitute(a), a2.Substitute(a), a3.Substitute(a), a4.Substitute(a), a6.Substitute(a));
    }

    public static EllGroup<GFelt> ToGF(this EllGroup<Rational> E, GFelt a) => E.ToGF(a.P).ToGF(a);

    public static EllGroup<ZnBigInt> ToZnBigInt(this EllGroup<Rational> E, BigInteger p)
    {
        var (a1, a2, a3, a4, a5) = E.Coefs;
        return new(a1.ToZnBigInt(p), a2.ToZnBigInt(p), a3.ToZnBigInt(p), a4.ToZnBigInt(p), a5.ToZnBigInt(p));
    }

    public static IEnumerable<EllPt<Rational>> SolveIntegralPoints(EllGroup<Rational> E, LogLevel lvl = LogLevel.Off)
    {
        var disc = E.Disc;
        var (A, B, C, _, _, _, _) = E.LongForm;
        var r = PrimesDec(BigInteger.Abs(disc.Num))
            .Aggregate(BigInteger.One, (acc, r) => acc * BigInteger.Pow(r.Key, r.Value / 2 + r.Value % 2));
        var divs = DividorsBigInt(16 * r).Where(y => (256 * disc.Num) % (y * y) == 0).Order()
            .Select(y => new Rational(y)).ToArray();

        var x = FG.QPoly();
        foreach (var y in divs.Prepend("0"))
        {
            var P = x.Pow(3) + A * x * x + B * x + C;
            var sols = IntFactorisation.FactorsQ(P - y.Pow(2));
            if (lvl == LogLevel.Level2)
                sols.Println($"Y = {y}, solve {y.Pow(2)} = {P}");

            var ellpts = sols.Where(e => e.Item1.Degree == 1).Select(e => new EllPt<Rational>(-e.Item1[0], y));
            foreach (var pt in ellpts)
            {
                yield return E.ConvertFromLong(pt);
                yield return E.ConvertFromLong(new(pt.X, -pt.Y));
            }
        }
    }

    public static void AddNewIntegralPoint(EllGroup<Rational> g, HashSet<EllPt<Rational>> set, EllPt<Rational> pt)
    {
        var sz = 0;
        while (set.Count != sz)
        {
            sz = set.Count;
            var tmp = set.Select(e => g.Op(e, pt)).Where(e => e.IsIntegral()).ToHashSet();
            set.UnionWith(tmp);
        }
    }

    public static (ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> pts, int[] abType)
        NagellLutzTorsionGroup(EllGroup<Rational> E, HashSet<EllPt<Rational>> ellpts, LogLevel lvl = LogLevel.Level1)
    {
        var set = new List<EllPt<Rational>>() { E.O };
        foreach (var pt in ellpts.Where(pt => E.ConvertToShort(pt).IsIntegral()))
        {
            var acc = pt;
            for (int i = 1; i <= 4; i++)
            {
                acc = E.Times(acc, 2);
                if (acc.IsO || !E.ConvertToShort(acc).IsIntegral())
                    break;
            }

            if (E.ConvertToShort(acc).IsIntegral())
                set.Add(pt);
        }

        var gEll = Group.Generate(E, set.ToArray());
        var abType = Group.AbelianGroupType(gEll);
        if (lvl != LogLevel.Off)
        {
            DisplayGroup.HeadElements(gEll);
            Console.WriteLine($"{gEll} TorsionGroup = {abType.Glue(" x ", "C{0}")}");
            Console.WriteLine();
        }

        var intPts = ellpts.Concat(gEll).Where(pt => pt.IsIntegral()).ToHashSet();
        return (gEll, intPts, abType);
    }

    public static (ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> pts, int[] abType)
        NagellLutzTorsionGroup(EllGroup<Rational> E, LogLevel lvl = LogLevel.Level1)
    {
        var ellpts = SolveIntegralPoints(E, lvl).ToHashSet();
        return NagellLutzTorsionGroup(E, ellpts, lvl);
    }

    static int Ord(int p, BigInteger A) => FactorMultiplicity(p, A).mul;

    static bool QuadRoot(Rational a0, Rational b0, Rational c0, int p)
    {
        var (a, b, c) = (a0.ToZnInt(p).K, b0.ToZnInt(p).K, c0.ToZnInt(p).K);
        if (a == 0)
            return b != 0 || c == 0;

        if (p == 2)
            return !(a == 1 && b == 1 && c == 1);

        var disc = AmodP(b * b - 4 * a * c, p);
        return disc == 0 || LegendreJacobi(disc, p) == 1;
    }

    static int NCubicRoots(Rational b0, Rational c0, Rational d0, int p)
    {
        var (b, c, d) = (b0.ToZnInt(p), c0.ToZnInt(p), d0.ToZnInt(p));
        var x = FG.ZPoly(p);
        var Q = x.Pow(3) + b * x * x + c * x + d;
        var a0 = NumberTheory.PrimitiveRootMod(p);
        var facts = IntFactorisation.FirrFsep(Q, x.KOne * a0);
        return facts.Count(f => f.g.Degree == 1);
    }

    public static TateAlgo TateAlgorithm(EllCoefs<Rational> E, int p)
    {
        var n = Ord(p, E.Disc.Num);
        if (n == 0)
            return new(p, n, "I0", 0, 1);

        var Etmp = E.Transform(0, 0, 0, 1);
        var (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;
        Rational r, s, t;
        int cp;
        if (p == 2)
        {
            if (b2.Mod(p).IsZero())
            {
                r = a4.Mod(p);
                t = (r * (1 + a2 + a4) + a6).Mod(p);
            }
            else
            {
                r = a3.Mod(p);
                t = (r + a4).Mod(p);
            }
        }
        else if (p == 3)
        {
            if (b2.Mod(p).IsZero())
                r = (-b6).Mod(p);
            else
                r = (-b2 * b4).Mod(p);

            t = (a1 * r + a3).Mod(p);
        }
        else
        {
            if (c4.Mod(p).IsZero())
                r = -b2 * new Rational(InvModPbezbigint(12, p));
            else
                r = -(c6 + b2 * c4) * new Rational(InvModPbezbigint(12 * c4.Num, p));

            t = -(a1 * r + a3) * (p + 1) / 2;
            r = r.Mod(p);
            t = t.Mod(p);
        }

        Etmp = Etmp.Transform(r, r.Zero, t, r.One);
        (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;

        if (!c4.Mod(p).IsZero())
        {
            if (QuadRoot(a1.One, a1, -a2, p))
                cp = n;
            else if (n % 2 == 0)
                cp = 2;
            else
                cp = 1;

            return new(p, n, $"In{n}", 1, cp);
        }

        if (!a6.Mod(p * p).IsZero())
            return new(p, n, "II", n, 1);

        if (!b8.Mod(p * p * p).IsZero())
            return new(p, n, "III", n - 1, 2);

        if (!b6.Mod(p * p * p).IsZero())
        {
            if (QuadRoot(a1.One, a3 / p, -a6 / (p * p), p))
                cp = 3;
            else
                cp = 1;

            return new(p, n, "IV", n - 2, cp);
        }

        if (p == 2)
        {
            s = a2.Mod(2);
            t = 2 * (a6 / 4).Mod(2);
        }
        else
        {
            s = -a1 * (p + 1) / 2;
            t = -a3 * (p + 1) / 2;
        }

        Etmp = Etmp.Transform(s.Zero, s, t, s.One);
        (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;

        var b = a2 / p;
        var c = a4 / (p * p);
        var d = a6 / p.Pow(3);
        var w = 27 * d * d - b * b * c * c + 4 * b.Pow(3) * d - 18 * b * c * d + 4 * c.Pow(3);
        var x = 3 * c - b * b;

        if (!w.Mod(p).IsZero())
            return new(p, n, "I*0", n - 4, 1 + NCubicRoots(b, c, d, p));
        else if (!x.Mod(p).IsZero())
        {
            if (p == 2)
                r = c;
            else if (p == 3)
                r = b * c;
            else
                r = (b * c - 9 * d) * new Rational(InvModPbezbigint(2 * x.Num, p));

            r = p * r.Mod(p);
            Etmp = Etmp.Transform(r, r.Zero, r.Zero, r.One);
            (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;

            var (m, mx, my) = (1, p * p, p * p);
            cp = 0;
            while (cp == 0)
            {
                var xa2 = a2 / p;
                var xa3 = a3 / my;
                var xa4 = a4 / (p * mx);
                var xa6 = a6 / (mx * my);
                if (!(xa3.Pow(2) + 4 * xa6).Mod(p).IsZero())
                {
                    if (QuadRoot(a1.One, xa3, -xa6, p))
                        cp = 4;
                    else
                        cp = 2;
                }
                else
                {
                    if (p == 2)
                        t = my * xa6;
                    else
                        t = my * (-xa3 * (p + 1) / 2).Mod(p);

                    Etmp = Etmp.Transform(t.Zero, t.Zero, t, t.One);
                    (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;

                    my = my * p;
                    m++;
                    xa2 = a2 / p;
                    xa3 = a3 / my;
                    xa4 = a4 / (p * mx);
                    xa6 = a6 / (mx * my);
                    if (!(xa4.Pow(2) - 4 * xa2 * xa6).Mod(p).IsZero())
                    {
                        if (QuadRoot(xa2, xa4, xa6, p))
                            cp = 4;
                        else
                            cp = 2;
                    }
                    else
                    {
                        if (p == 2)
                            r = mx * (xa6 * xa2).Mod(2);
                        else
                            r = mx * (-xa4 * new Rational(InvModPbezbigint(2 * xa2.Num, p))).Mod(p);

                        Etmp = Etmp.Transform(r, r.Zero, r.Zero, r.One);
                        (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;
                        mx *= p;
                        m++;
                    }
                }
            }

            return new(p, n, $"I*n{m}", n - m - 4, cp);
        }
        else
        {
            var rp = r.Zero;
            if (p == 3)
                rp = -d;
            else
                rp = -b * new Rational(InvModPbezbigint(3, p));

            r = p * rp.Mod(p);
            Etmp = Etmp.Transform(r, r.Zero, r.Zero, r.One);
            (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;
            var x3 = a3 / p.Pow(2);
            var x6 = a6 / p.Pow(4);
            if (!(x3.Pow(2) + 4 * x6).Mod(p).IsZero())
            {
                if (QuadRoot(a1.One, x3, -x6, p))
                    cp = 3;
                else
                    cp = 1;

                return new(p, n, "IV*", n - 6, cp);
            }
            else
            {
                if (p == 2)
                    t = x6;
                else
                    t = x3 * (p + 1) / 2;

                t = -p * p * t.Mod(p);
                Etmp = Etmp.Transform(t.Zero, t.Zero, t, t.One);
                (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, _) = Etmp.ModelAndInvariants;

                if (!a4.Mod(p.Pow(4)).IsZero())
                    return new(p, n, "III*", n - 7, 2);
                else if (!a6.Mod(p.Pow(6)).IsZero())
                    return new(p, n, "II*", n - 8, 1);
                else
                    return TateAlgorithm(Etmp.Transform(0, 0, 0, p), p);
            }
        }
    }

    public static int EllAp(EllGroup<Rational> E, int p)
    {
        if (p <= 3)
        {
            var (_a1, _a2, _a3, _a4, _a6) = E.Coefs;
            var (a1, a2, a3, a4, a6) = (_a1.ToZnInt(p), _a2.ToZnInt(p), _a3.ToZnInt(p), _a4.ToZnInt(p), _a6.ToZnInt(p));
            var card = 1 + p.Range().Select(k => new ZnInt(p, k)).ToArray().Grid2D().Select(e => (x: e.t1, y: e.t2))
                .Count(
                    e => (e.y.Pow(2) + a1 * e.x * e.y + a3 * e.y).Equals(e.x.Pow(3) + a2 * e.x * e.x + a4 * e.x + a6));

            return p + 1 - card;
        }

        var (A, B, _, _, _, _) = E.ShortForm;
        var (a, b) = (A.ToZnBigInt(p).Unsigned, B.ToZnBigInt(p).Unsigned);
        return -p.Range().Select(x => (int)LegendreJacobiBigint((BigInteger.Pow(x, 3) + a * x + b) % p, p))
            .Sum(k => k <= 1 ? k : -1);
    }

    static List<BigInteger> CandidatsN(double pmin, double pmax, BigInteger L)
    {
        var n0 = (int)(pmin / (double)L);
        var n1 = (int)(pmax / (double)L) + 1;
        return (n1 - n0 + 1).SeqLazy(n0).Where(e => e > 0 && L % e == 0)
            .Select(e => e * L).Where(n => (double)n > pmin && (double)n < pmax).ToList();
    }

    static EllPt<ZnBigInt> RandPt(EllGroup<ZnBigInt> Ep, HashSet<BigInteger> set)
    {
        var p = Ep.ShortForm.A.Mod;
        while (true)
        {
            var xP = new ZnBigInt(p, DistributionExt.Dice(BigInteger.Zero, p));
            var y2P = xP.Pow(3) + Ep.ShortForm.A * xP + Ep.ShortForm.B;
            if (LegendreJacobiBigint(y2P.Unsigned, p) != 1)
                continue;

            var yP = xP.One * NumberTheory.SqrtModANTV1(y2P.K, p);
            if (set.Add(xP.Unsigned))
                return Ep.ConvertFromShort(new(xP, yP));
        }
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

        var g = 1000.SeqLazy().Select(_ => DistributionExt.Dice(BigInteger.One * 2, p - 1))
            .Where(g => LegendreJacobiBigint(g, p) != 1)
            .Select(g => new ZnBigInt(p, g))
            .Where(g => !(4 * (E.ShortForm.A.ToZnBigInt(p) * g.Pow(2)).Pow(3) +
                          27 * (E.ShortForm.B.ToZnBigInt(p) * g.Pow(3)).Pow(2)).IsZero())
            .First();
        var Ep = E.ToZnBigInt(p);
        // Mestre's Theorem, quadratic twist of Ep
        var Ep_ = new EllGroup<ZnBigInt>(Ep.ShortForm.A * g.Pow(2), Ep.ShortForm.B * g.Pow(3));
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
                Console.WriteLine($"Ep = {Ep.ShortFormStr} Ap = {Ap} Case = [{caseNum}]#{outputNum}");

            return Ap;
        }
        else if (outputNum == 2)
        {
            var Ap = (int)cands_.Where(NL => (p - 1) % (NL / L_) == 0).Select(NL => -(p + 1 - NL)).First();
            if (log != LogLevel.Off)
                Console.WriteLine($"Ep = {Ep.ShortFormStr} Ap = {Ap} Case = [{caseNum}]#2");

            return Ap;
        }
        else
            throw new();
    }

    public static Dictionary<int, int> EllAn(EllGroup<Rational> E, Rational N)
    {
        var nmax = 2 * Rational.Sqrt(N).Num;
        var ellAn = Primes10000.Where(p => p <= nmax).ToDictionary(p => p, p => EllAp(E, p));
        ellAn[0] = 0;
        ellAn[1] = 1;
        foreach (var p in ellAn.Keys.Where(p => p > 1).Order().ToArray())
        {
            for (int i = 2; p.Pow(i) <= nmax; i++)
            {
                var t = N.Mod(p).IsZero() ? 0 : 1;
                var e1 = ellAn[p];
                var e2 = ellAn[p.Pow(i - 1)];
                var e3 = ellAn[p.Pow(i - 2)];
                var e4 = e2 * e1 - t * p * e3;
                ellAn[p.Pow(i)] = e4;
            }
        }

        for (int i = 0; i <= nmax; i++)
        {
            if (ellAn.ContainsKey(i))
                continue;

            var dec = PrimesDec(i);
            ellAn[i] = dec.Select(e => ellAn[e.Key.Pow(e.Value)]).Aggregate(1, (acc, aj) => acc * aj);
        }

        return ellAn;
    }

    public static double factorial(int k)
    {
        var f = 1.0;
        for (int i = 1; i <= k; i++)
            f *= i;

        return f;
    }

    static double[] PrCoefs =
    [
        1, -0.57721566490153286060651209008240243104, 0.98905599532797255539539565150063470794,
        -0.90747907608088628901656016735627511493, 0.98172808683440018733638029402185085036,
        -0.98199506890314520210470141379137467551, 0.99314911462127619315386725332865849803,
        -0.99600176044243153397007841966456668673, 0.99810569378312892197857540308836723751,
        -0.99902526762195486779467805964888808852, 0.99951565607277744106705087759437019442,
        -0.99975659750860128702584244914060923598, 0.99987827131513327572617164259000321937,
        -0.99993906420644431683585223136895513183, 0.99996951776348210449861140509195350724,
        -0.99998475269937704874370963172444753831, 0.99999237447907321585539509450510782580,
        -0.99999618658947331202896495779561431377, 0.99999809308113089205186619151459489769,
        -0.99999904646891115771748687947054372628, 0.99999952321060573957523929299106456813,
        -0.99999976159734438057092470106258744744, 0.99999988079601916841665041840424924048,
        -0.99999994039712498374586288797675081780, 0.99999997019826758235557449619251141976
    ];

    private const double EPS = 0.0001;

    static double Pr(int r, double x) => (r + 1).SeqLazy().Sum(k => PrCoefs[r - k] * double.Pow(x, k) / factorial(k));

    public static double G(int r, double x)
    {
        if (x > 30)
            return 0;

        var pr = Pr(r, double.Log(1 / x));
        var tmp = pr - 1;
        var (n, sgn) = (1, -(-1).Pow(r));
        while (double.Abs(pr - tmp) > 1e-12)
        {
            tmp = pr;
            pr += sgn * double.Pow(x, n) / (double.Pow(n, r) * factorial(n));
            ++n;
            sgn *= -1;
        }

        return 2 * pr;
    }

    public static double Lr(int r, double X, Dictionary<int, int> ellAn) =>
        factorial(r) * ellAn.Keys.Max().SeqLazy(1).Sum(n => ellAn[n] * G(r, X * n) / n);

    public static double L0Approx(double A, double X, Dictionary<int, int> ellAn) =>
        ellAn.Keys.Max().SeqLazy(1).Sum(n => ellAn[n] * (G(0, A * X * n) - G(0, X * n / A)) / n);

    public static (EllGroup<Rational> Ell, Rational N, TateAlgo[] tate) EllTateAlgorithm(EllCoefs<Rational> E)
    {
        var Ell = E.ToEllGroup();
        var dec = PrimesDec(E.Disc.Absolute.Num);
        var tate = dec.Keys.Select(p => TateAlgorithm(E, p)).ToArray();
        var N = tate.Select(e => new Rational(e.p).Pow(e.fp)).Aggregate((pi, pj) => pi * pj);
        return (Ell, N, tate);
    }

    public static (int rank, double Lr, Rational N, Dictionary<int, int> ellAn, TateAlgo[], EllGroup<Rational> Ell)
        AnalyticRank(EllCoefs<Rational> E)
    {
        var (Ell, N, tate) = EllTateAlgorithm(E);
        var ellAn = EllAn(Ell, N);

        var X = 2 * double.Pi / double.Sqrt(N);

        var rank = 0;
        var L0 = Lr(0, X, ellAn);
        var L1 = Lr(1, X, ellAn);
        double L = 0.0;

        if (double.Abs(L0) > EPS && double.Abs(L1) > EPS)
        {
            var L0a = L0Approx(1.1, X, ellAn);
            if (double.Abs(L0a) < EPS)
            {
                rank = 1;
                L = L1;
            }
            else
                L = L0;
        }
        else
        {
            if (double.Abs(L0) < EPS)
                rank = 2;
            else
                rank = 3;

            for (; rank < 8; rank += 2)
            {
                var lr = Lr(rank, X, ellAn);
                L = lr;
                if (double.Abs(lr) > EPS)
                    break;
            }
        }

        return (rank, L, N, ellAn, tate, Ell);
    }


    public static (int rank, double Lr, Rational N, Dictionary<int, int> ellAn, TateAlgo[] tate, EllGroup<Rational> Ell)
        AnalyticRank(BigInteger[] curve) => AnalyticRank(EllCoefs(curve));

    public static KPoly<KPoly<K>> Simplify<K>(KPoly<KPoly<K>> P, KPoly<KPoly<K>> F, KPoly<KPoly<K>> G)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (P.Degree < F.Degree)
            return P;

        var (quo, rem) = P.Div(F);
        return Simplify(quo * G + rem, F, G);
    }

    public static (EllPt<FracPoly<FracPoly<K>>> P, EllPt<FracPoly<FracPoly<K>>> nP)
        NPt<K>(int n, Dictionary<int, KPoly<KPoly<K>>> psi, KPoly<KPoly<K>> R0, KPoly<KPoly<K>> R1)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var Y = FG.KFracPoly('Y', FG.KFracPoly(R1.KOne.X).X).X;
        var X = Y.KOne.X + Y.Zero;

        var tmp1 = Simplify(psi[n - 1] * psi[n + 1], R0, R1).ToFrac(Y) /
                   Simplify(psi[n].Pow(2), R0, R1).ToFrac(Y);
        var nPX = X - tmp1;

        var denom = n % 2 == 0
            ? Y * Simplify(4 * psi[n].Pow(3), R0, R1).ToFrac(Y)
            : Simplify(4 * R1.X * psi[n].Pow(3), R0, R1).ToFrac(Y);
        var tmp2 = Simplify(psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2), R0, R1).ToFrac(Y);
        var nPY = tmp2 / denom;

        return (new EllPt<FracPoly<FracPoly<K>>>(X, Y), new EllPt<FracPoly<FracPoly<K>>>(nPX, nPY));
    }

    public static (KPoly<KPoly<K>> R0, KPoly<KPoly<K>> R1, Dictionary<int, KPoly<KPoly<K>>> psi,
        Dictionary<int, KPoly<K>> f)
        DivisionPolynomial<K>(EllGroup<K> E, int nmax) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = FG.KPoly('X', E.a1);
        var Y = FG.KPoly('Y', x);
        var X = x * Y.One;
        var (X2, X3, X4, X5, X6) = 5.SeqLazy(2).Select(i => X.Pow(i)).Deconstruct();
        var (A, B) = (E.ShortForm.A * x.One * Y.One, E.ShortForm.B * x.One * Y.One);
        var eCoefs = new EllCoefs<K>(E.a1, E.a2, E.a3, E.a4, E.a5);
        var (a1, a2, a3, a4, a6) = (E.a1 * x.One * Y.One, E.a2 * x.One * Y.One, E.a3 * x.One * Y.One,
            E.a4 * x.One * Y.One,
            E.a5 * x.One * Y.One);
        var (b2, b4, b6, b8) = (eCoefs.b2 * x.One * Y.One, eCoefs.b4 * x.One * Y.One, eCoefs.b6 * x.One * Y.One,
            eCoefs.b8 * x.One * Y.One);

        var R0 = Y.Pow(2) + a1 * X * Y + a3 * Y;
        var R1 = X3 + a2 * X2 + a4 * X + a6;

        var psi = new Dictionary<int, KPoly<KPoly<K>>>();
        (psi[0], psi[1], psi[2]) = (Y.Zero, Y.One, 2 * Y + a1 * X + a3);
        psi[3] = 3 * X4 + b2 * X3 + 3 * b4 * X2 + 3 * b6 * X + b8;
        psi[4] = psi[2] * (2 * X6 + b2 * X5 + 5 * b4 * X4 + 10 * b6 * X3 + 10 * b8 * X2 + (b2 * b8 - b4 * b6) * X +
                           (b4 * b8 - b6 * b6));

        for (int n = 2; n < nmax / 2; n++)
        {
            psi[2 * n] = Simplify(psi[n] * (psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2)) / psi[2],
                R0, R1);
            psi[2 * n + 1] = Simplify(psi[n + 2] * psi[n].Pow(3) - psi[n + 1].Pow(3) * psi[n - 1], R0, R1);
        }

        var f = psi.ToDictionary(e => e.Key, e => e.Key % 2 == 0 ? (e.Value / Y)[0].Clone() : e.Value[0].Clone());
        f[2] = NPt(2, psi, R0, R1).nP.X.Num[0].Denom;
        return (R0, R1, psi, f);
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
        var facts = IntFactorisation.MusserSFF(P)
            .SelectMany(f => IntFactorisation.CantorZassenhausAECF(f.g, a, q))
            .ToArray();

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
            var sqrts0 = seq0.Select(e => (e.x, e.b, e.c, e.delta, sqrtDelta: NumberTheory.SqrtModANTV1(e.delta, g0)));
            return (sqrts0.SelectMany(e => new[]
                {
                    new EllPt<GFelt>(e.x, (-e.b + e.sqrtDelta) / 2),
                    new EllPt<GFelt>(e.x, (-e.b - e.sqrtDelta) / 2)
                })
                .ToHashSet(), g0);
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
            var (X, a0) = FG.EPolyXc(r, delta.Poly.x);
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
            var sqrts1 = seq1.Select(e => (e.x, e.b, e.c, e.delta, sqrtDelta: NumberTheory.SqrtModANTV1(e.delta, g1)));
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
        var Efq1 = Efq.ToGF(g.X);
        if (pts.Any(pt => !Efq1.Contains(pt.X, pt.Y)))
            throw new("#2");

        var orders = pts.ToDictionary(pt => pt, pt => Group.Cycle(Efq1, pt));
        var e1 = orders.First(e => e.Value.Count == l).Key;
        var e2 = orders.First(e => e.Value.Count == l && !e.Value.Keys.Contains(e1)).Key;
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
}