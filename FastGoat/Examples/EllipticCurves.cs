using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class EllipticCurves
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

    static int[][] AbSubTypes(int[] type)
    {
        var all = type.Select(t => IntExt.DividorsInt(t).Order().ToArray()).MultiLoop()
            .Select(l => l.Order().ToArray())
            .ToHashSet(new SequenceEquality<int>())
            .Select(l => l.Where(e => e != 1).ToArray())
            .Where(l => l.Length != 0)
            .Append([1])
            .OrderBy(l => l.Length)
            .ThenBy(l => l, Comparer<int[]>.Create((l0, l1) => l0.SequenceCompareTo(l1)))
            .ToArray();

        return all;
    }

    static T[] CurveArray<T>(params T[] ts) where T : struct, IElt<T>, IRingElt<T>
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

    static ConcreteGroup<EllPt<ZnBigInt>> EllFp(EllGroup<ZnBigInt> E)
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

    static int[] EllFpType(EllGroup<ZnBigInt> E, LogLevel lvl = LogLevel.Off)
    {
        if (E.Disc.IsZero())
            return [];

        var gEll = EllFp(E);
        if (lvl == LogLevel.Level2)
            DisplayGroup.HeadElements(gEll);

        var abType = Group.AbelianGroupType(gEll);
        Console.WriteLine($"{gEll} ~ {abType.Glue(" x ", "C{0}")}");
        if (gEll.Any(e => !e.IsO && !E.Contains(e.X, e.Y)))
            throw new();

        return abType;
    }

    static EllGroup<ZnBigInt> EllRational2ZnBigInt(EllGroup<Rational> E, int p)
    {
        var (a1, a2, a3, a4, a5) = E.Coefs;
        return new(a1.ToZnBigInt(p), a2.ToZnBigInt(p), a3.ToZnBigInt(p), a4.ToZnBigInt(p), a5.ToZnBigInt(p));
    }

    static void EllTors(EllGroup<Rational> E, int nbPrimes = 10, LogLevel lvl = LogLevel.Off)
    {
        Console.WriteLine($"#### Start {E}");
        var disc = 16 * E.Disc;
        var allTypes = IntExt.Primes10000.Where(p => p > 3 && (2 * disc) % p != 0).Take(nbPrimes)
            .Select(p => EllRational2ZnBigInt(E, p))
            .Select(Ep => EllFpType(Ep, lvl))
            .ToHashSet(new SequenceEquality<int>())
            .ToArray();

        var allSubTypes = allTypes.Select(l => AbSubTypes(l.ToArray()).ToHashSet(new SequenceEquality<int>()))
            .ToArray();
        var set = new HashSet<IEnumerable<int>>(new SequenceEquality<int>());
        foreach (var sub in allSubTypes)
        {
            if (set.Count == 0)
            {
                set = sub;
                continue;
            }

            set.IntersectWith(sub);
        }

        allTypes.Println(e => e.Glue(" x ", "C{0}"), $"Morphisms {E} ->");
        set.Select(l => l.ToArray()).Println(e => e.Glue(" x ", "C{0}"), "Intersections subgroups");
        var tor = set.MaxBy(l => l.Aggregate(1, (acc, i) => acc * i));
        Console.WriteLine($"{E} Torsion ~-> {tor!.Descending().Glue(" x ", "C{0}")}");
        Console.WriteLine();
    }

    static void EllTors(BigInteger a, BigInteger b, int nbPrimes = 10, LogLevel lvl = LogLevel.Off)
    {
        EllTors(new(new Rational(a), new Rational(b)), nbPrimes, lvl);
    }

    static IEnumerable<EllPt<Rational>> SolveIntegralPoints(EllGroup<Rational> E, LogLevel lvl = LogLevel.Off)
    {
        var disc = E.Disc;
        var (A, B, C, _, _) = E.LongForm;
        // Console.WriteLine($"{E} ~ {E.LongFormStr} ~ {E.ShortFormStr}");
        // var (a1, a2, a3, a4, a5) = E.Coefs;
        // var (X, Y) = Ring.Polynomial(Rational.KOne(), "X", "Y").Deconstruct();
        // var F = X.Pow(3) + a2 * X.Pow(2) + a4 * X + a5 - (Y.Pow(2) + a1 * X * Y + a3 * Y);
        // var disc1 = Ring.Discriminant(Ring.Discriminant(F, Y) / 16, X).ConstTerm / 16;
        // Console.WriteLine(new { disc, disc1, div = disc / disc1 });
        var r = IntExt.PrimesDec(BigInteger.Abs(disc.Num))
            .Aggregate(BigInteger.One, (acc, r) => acc * BigInteger.Pow(r.Key, r.Value / 2 + r.Value % 2));
        var divs = IntExt.DividorsBigInt(16 * r).Where(y => (256 * disc.Num) % (y * y) == 0).Order()
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

    static void AddNewIntegralPoint(EllGroup<Rational> g, HashSet<EllPt<Rational>> set, EllPt<Rational> pt)
    {
        var sz = 0;
        while (set.Count != sz)
        {
            sz = set.Count;
            var tmp = set.Select(e => g.Op(e, pt)).Where(e => e.IsIntegral()).ToHashSet();
            set.UnionWith(tmp);
        }
    }

    static (ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> pts, int[] abType)
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

    static (ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> pts, int[] abType)
        NagellLutzTorsionGroup(EllGroup<Rational> E, LogLevel lvl = LogLevel.Level1)
    {
        var ellpts = SolveIntegralPoints(E, lvl).ToHashSet();
        return NagellLutzTorsionGroup(E, ellpts, lvl);
    }

    static (ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> pts, int[] abType)
        NagellLutzTorsionGroup(BigInteger[] curve, LogLevel lvl = LogLevel.Level1)
    {
        var (a1, a2, a3, a4, a5) = CurveArray(curve.Select(e => new Rational(e)).ToArray()).Deconstruct();
        return NagellLutzTorsionGroup(new EllGroup<Rational>(a1, a2, a3, a4, a5), lvl);
    }

    static (ConcreteGroup<EllPt<Rational>> gEll, HashSet<EllPt<Rational>> pts, int[] abType)
        NagellLutzTorsionGroup(BigInteger[] curve, (BigInteger x, BigInteger y)[] ellPts,
            LogLevel lvl = LogLevel.Level1)
    {
        var (a1, a2, a3, a4, a5) = CurveArray(curve.Select(e => new Rational(e)).ToArray()).Deconstruct();
        var pts = ellPts.Select(e => new EllPt<Rational>($"{e.x}", $"{e.y}")).ToHashSet();
        return NagellLutzTorsionGroup(new EllGroup<Rational>(a1, a2, a3, a4, a5), pts, lvl);
    }

    static void Transform((Polynomial<Rational, Xi> lhs, Polynomial<Rational, Xi> rhs) e,
        TorsionMeth meth = TorsionMeth.Both)
    {
        var ((y, _), (x, X)) = (-e.lhs + e.rhs).IndeterminatesAndVariables.Deconstruct();
        var ind = X.Indeterminates;
        var (xm, ym) = (new Monom<Xi>(ind, x), new Monom<Xi>(ind, y));
        var (xym, x2m) = (xm.Mul(ym), xm.Pow(2));
        var (a1, a2, a3, a4, a5) = (e.lhs[xym], e.rhs[x2m], e.lhs[ym], e.rhs[xm], e.rhs.ConstTerm);
        var E = new EllGroup<Rational>(a1, a2, a3, a4, a5);
        var P1 = X.Pow(3) + E.LongForm.A * X.Pow(2) + E.LongForm.B * X + E.LongForm.C;
        var P2 = X.Pow(3) + E.ShortForm.A * X + E.ShortForm.B;
        Console.WriteLine($"Elliptic curve      {e.lhs} = {e.rhs}");
        Console.WriteLine($"Simplified form     y^2 = {P1}");
        Console.WriteLine($"Simplified form     y^2 = {P2}");

        if ((meth & TorsionMeth.Fp) == TorsionMeth.Fp)
            EllTors(E);
        if ((meth & TorsionMeth.NagellLutz) == TorsionMeth.NagellLutz)
            NagellLutzTorsionGroup(E);

        Console.WriteLine();
    }

    static (int y, bool sol) ApproxSolver(int x, int n)
    {
        var y2 = BigInteger.Pow(x, 3) + n;
        if (y2 < 0)
            return (0, false);

        var ya = BigInteger.Parse($"{double.Ceiling(double.Sqrt((double)y2))}");
        while (ya * ya < y2 + 1)
        {
            if (ya * ya == y2)
                return ((int)ya, true);

            ya++;
        }

        return (0, false);
    }

    public static int SchoofEllPtsCount(BigInteger a, BigInteger b, int p)
    {
        return p + 1 + p.Range()
            .Select(x => IntExt.LegendreJacobiBigint((BigInteger.ModPow(x, 3, p) + a * x + b) % p, p))
            .Sum(k => k <= 1 ? (int)k : -1);
    }

    public static int EllRank(BigInteger a, BigInteger b, int n = 500, bool show = true)
    {
        var r = 1.0;
        var (sumX, sumY, sumX2, sumXY) = (0.0, 0.0, 0.0, 0.0);
        foreach (var p in IntExt.Primes10000.Take(n))
        {
            r *= 1.0 * SchoofEllPtsCount(a, b, p) / p;

            var (x, y) = (double.Log(double.Log(p)), double.Log(r));
            sumX += x;
            sumY += y;
            sumX2 += x * x;
            sumXY += x * y;
        }

        var A = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
        var B = (sumY * sumX2 - sumX * sumXY) / (n * sumX2 - sumX * sumX);
        var rk = (int)double.Round(A);
        if (show)
        {
            Console.WriteLine($"Regr Line y = {A:f4} * x + {B:f4}");
            Console.WriteLine($"Rank(E[{a},{b}](Q)) = {rk}");
            Console.WriteLine();
        }

        return rk;
    }

    static void symbWeierstrassForm()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.GrLex,
            "a1", "a2", "a3", "a4", "a5", "d1", "d2", "A", "B", "C", "x", "y");
        var (a1, a2, a3, a4, a5, d1, d2) = Ring.EPolynomial(xis.Take(7).ToArray()).Deconstruct();
        var (A, B, C, X, Y) = Ring.EPolynomial(xis.Skip(7).ToArray()).Deconstruct();
        var (x, y) = (X.Num.ExtractIndeterminate, Y.Num.ExtractIndeterminate);

        var eqEll = (X.Pow(3) + a2 * X.Pow(2) + a4 * X + a5) - (Y.Pow(2) + a1 * X * Y + a3 * Y);

        var eqLong1 = eqEll.Substitute(Y - (a1 * X + a3) / 2, y);
        Console.WriteLine(eqEll);
        Console.WriteLine(eqLong1);
        Ring.Decompose(eqLong1.Num, x).Item1.Println();
        // (-a1*x*y + a2*x^2 + x^3 - a3*y + a4*x - y^2 + a5)
        // (1/4*a1^2*x^2 + 1/2*a1*a3*x + a2*x^2 + x^3 - 1/4*a2^2 + 1/2*a2*a3 + a2*y - a3*y + a4*x - y^2 + a5)
        // Lines
        //     [1, -1/4*a2^2 + 1/2*a2*a3 + a2*y - a3*y - y^2 + a5]
        //     [x, 1/2*a1*a3 + a4]
        //     [x^2, 1/4*a1^2 + a2]
        //     [x^3, 1]
        // 

        var eqLong2 = (X.Pow(3) + A * X * X + B * X + C) - Y.Pow(2);
        var eqLong3 = (d1.Pow(3) * eqLong2.Substitute(X / d1, x)).Substitute(Y / d2, y);
        Console.WriteLine(eqLong2);
        Console.WriteLine(eqLong3);
        Console.WriteLine(eqLong3.Num);
        Console.WriteLine(eqLong3.Denom);
        Ring.Decompose(eqLong3.Num, x).Item1.Println();
        Console.WriteLine();
        // (A*x^2 + x^3 + B*x - y^2 + C)
        // (-C*d1^3*d2^2 - B*d1^2*d2^2*x - A*d1*d2^2*x^2 + d1^3*y^2 - d2^2*x^3)/(-d2^2)
        // -C*d1^3*d2^2 - B*d1^2*d2^2*x - A*d1*d2^2*x^2 + d1^3*y^2 - d2^2*x^3
        // -d2^2
        // Lines
        //     [1, -C*d1^3*d2^2 + d1^3*y^2]
        //     [x, -B*d1^2*d2^2]
        //     [x^2, -A*d1*d2^2]
        //     [x^3, -d2^2]
        // 

        var eqShort1 = eqEll.Substitute((Y - (a1 * X + a3) / 2, y), (X - a1.Pow(2) / 12 - a2 / 3, x));
        Console.WriteLine(eqEll);
        Console.WriteLine(eqShort1);
        Ring.Decompose(eqShort1.Num, x).Item1.Println();
        // (-a1*x*y + a2*x^2 + x^3 - a3*y + a4*x - y^2 + a5)
        // (1/864*a1^6 + 1/72*a1^4*a2 - 1/48*a1^4*x - 1/24*a1^3*a3 + 1/18*a1^2*a2^2 - 1/6*a1^2*a2*x - 1/12*a1^2*a4 - 1/6*a1*a2*a3 + 1/2*a1*a3*x + 2/27*a2^3 - 1/3*a2^2*x + x^3 - 1/3*a2*a4 + 1/4*a3^2 + a4*x - y^2 + a5)
        // Lines
        //     [1, 1/864*a1^6 + 1/72*a1^4*a2 - 1/24*a1^3*a3 + 1/18*a1^2*a2^2 - 1/12*a1^2*a4 - 1/6*a1*a2*a3 + 2/27*a2^3 - 1/3*a2*a4 + 1/4*a3^2 - y^2 + a5]
        //     [x, -1/48*a1^4 - 1/6*a1^2*a2 + 1/2*a1*a3 - 1/3*a2^2 + a4]
        //     [x^2, 0]
        //     [x^3, 1]
        // 

        var eqShort2 = (X.Pow(3) + A * X + B) - Y.Pow(2);
        var eqShort3 = (d1.Pow(3) * eqShort2.Substitute(X / d1, x)).Substitute(Y / d2, y);
        Console.WriteLine(eqShort2);
        Console.WriteLine(eqShort3);
        Console.WriteLine(eqShort3.Num);
        Console.WriteLine(eqShort3.Denom);
        Ring.Decompose(eqShort3.Num, x).Item1.Println();
        // (-B*d1^3*d2^2 - A*d1^2*d2^2*x + d1^3*y^2 - d2^2*x^3)/(-d2^2)
        // -B*d1^3*d2^2 - A*d1^2*d2^2*x + d1^3*y^2 - d2^2*x^3
        // -d2^2
        // Lines
        //     [1, -B*d1^3*d2^2 + d1^3*y^2]
        //     [x, -A*d1^2*d2^2]
        //     [x^2, 0]
        //     [x^3, -d2^2]
        // 

        Console.WriteLine();
    }

    // Elliptic Curve Discrete Logarithm Problem
    static void ECDLP(int p, (int x, int y) P0, (int x, int y) Q0, int[] curve)
    {
        var (a1, a2, a3, a4, a5) = CurveArray(curve.Select(i => new ZnInt(p, i)).ToArray()).Deconstruct();
        var E = new EllGroup<ZnInt>(a1, a2, a3, a4, a5);
        var nb = SchoofEllPtsCount(E.ShortForm.A.K, E.ShortForm.B.K, p);
        Console.WriteLine($"|{E}| = {nb}");

        var P = new EllPt<ZnInt>(new(p, P0.x), new(p, P0.y));
        var Q = new EllPt<ZnInt>(new(p, Q0.x), new(p, Q0.y));

        GlobalStopWatch.AddLap();
        var k = Group.BSGS(E, P, Q, nb);
        var Q1 = E.Times(P, k);
        Console.WriteLine($"P={P} Q={Q1} {k}xP=Q");
        GlobalStopWatch.Show();
        Console.WriteLine();
        if (!Q.Equals(Q1))
            throw new();
    }

    public static void Example1()
    {
        var E = new EllGroup<Rational>("-36", "0");
        var O = E.O;
        EllPt<Rational> P = ("-3", "9");
        EllPt<Rational> Q = ("-2", "8");
        Console.WriteLine(new { O, P, Q });
        Console.WriteLine($"-P = {E.Invert(P)}");
        Console.WriteLine($"P + Q = {E.Op(P, Q)}");
        Console.WriteLine($"2P = {E.Times(P, 2)}");
        Console.WriteLine($"2Q = {E.Times(Q, 2)}");
        Console.WriteLine($"2P + 2Q = {E.Op(E.Times(P, 2), E.Times(Q, 2))}");
        Console.WriteLine($"2(P + Q) = {E.Times(E.Op(P, Q), 2)}");
    }

    // [-36,0]
    // https://www.lmfdb.org/EllipticCurve/Q/576/c/3
    // 
    // [1,0]
    // https://www.lmfdb.org/EllipticCurve/Q/64/a/4
    // 
    // [0,3]
    // https://www.lmfdb.org/EllipticCurve/Q/3888/i/2
    // 
    // [-43,166]
    // https://www.lmfdb.org/EllipticCurve/Q/26/b/2
    public static void Example2Fp()
    {
        EllTors(-36, 0, lvl: LogLevel.Level2); // Ell[-36,0](Q) Torsion = C2 x C2
        EllTors(0, 3); // Ell[0,3](Q) Torsion = C1
        EllTors(1, 0); // Ell[1,0](Q) Torsion = C2
        EllTors(0, 1); // Ell[0,1](Q) Torsion = C6
        EllTors(-43, 166, nbPrimes: 20); // Ell[-43,166](Q) Torsion = C7
    }

    public static void Example3NagellLutz()
    {
        NagellLutzTorsionGroup([-36, 0], lvl: LogLevel.Level2); // Ell[-36,0](Q) Torsion = C2 x C2
        NagellLutzTorsionGroup([0, 3]); // Ell[0,3](Q) Torsion = C1
        NagellLutzTorsionGroup([1, 0]); // Ell[1,0](Q) Torsion = C2
        NagellLutzTorsionGroup([0, 1]); // Ell[0,1](Q) Torsion = C6
        NagellLutzTorsionGroup([-43, 166]); // Ell[-43,166](Q) Torsion = C7
    }
    // |Ell[-43,166](Q)| = 7
    // Type        AbelianGroup
    // BaseGroup   Ell[-43,166](Q)
    // 
    // Elements
    // (1)[1] = O
    // (2)[7] = (-5,-16)
    // (3)[7] = (-5,16)
    // (4)[7] = (3,-8)
    // (5)[7] = (3,8)
    // (6)[7] = (11,-32)
    // (7)[7] = (11,32)
    // 
    // Ell[-43,166](Q) Torsion = C7

    public static void Example4TransformCurve()
    {
        Logger.SetOff();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        symbWeierstrassForm();

        var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();

        // – E1 : y^2 + 7xy = x^3 + 16x ;
        var e1 = (y.Pow(2) + 7 * x * y, x.Pow(3) + 16 * x);
        // – E2 : y^2 + xy − 5y = x^3 − 5x^2 ;
        var e2 = (y.Pow(2) + x * y - 5 * y, x.Pow(3) - 5 * x.Pow(2));
        // – E3 : y^2 − y = x^3 − x^2 ;
        var e3 = (y.Pow(2) - y, x.Pow(3) - x.Pow(2));
        // – E4 : y^2 + xy + y = x^3 − x^2 − 14x + 29 ;
        var e4 = (y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 14 * x + 29);
        // – E5 : y^2 + xy = x^3 − 45x + 81 ;
        var e5 = (y.Pow(2) + x * y, x.Pow(3) - 45 * x + 81);
        // – E6 : y^2 + 43xy − 210y = x^3 − 210x^2 ;
        var e6 = (y.Pow(2) + 43 * x * y - 210 * y, x.Pow(3) - 210 * x.Pow(2));
        // – E7 : y^2 + 5xy − 6y = x^3 − 3x^2 ;
        var e7 = (y.Pow(2) + 5 * x * y - 6 * y, x.Pow(3) - 3 * x.Pow(2));
        // – E8 : y^2 + 17xy − 120y = x^3 − 60x^2 ;
        var e8 = (y.Pow(2) + 17 * x * y - 120 * y, x.Pow(3) - 60 * x.Pow(2));
        // – E9 : y^2 + xy = x^3 − 1070x + 7812 ;
        var e9 = (y.Pow(2) + x * y, x.Pow(3) - 1070 * x + 7812);

        Transform(e1);
        Transform(e2);
        Transform(e3);
        Transform(e4);
        Transform(e5);
        Transform(e6);
        Transform(e7);
        Transform(e8);
        Transform(e9);
    }

    public static void Example5FromLMFDB()
    {
        Logger.SetOff();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();
        GlobalStopWatch.Restart();

        // Torsion C1
        Transform((y.Pow(2), x.Pow(3) - 4 * x - 4), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + y, x.Pow(3) + x), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 29 * x - 53), TorsionMeth.NagellLutz);

        // Torsion C2
        Transform((y.Pow(2), x.Pow(3) - 11 * x - 14), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) + x.Pow(2) + x), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - 14 * x - 64), TorsionMeth.NagellLutz);

        // Torsion C3
        Transform((y.Pow(2), x.Pow(3) + 4), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + y, x.Pow(3) + x.Pow(2) + x - 1), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) + 6 * x - 28), TorsionMeth.NagellLutz);

        // Torsion C4, C2 x C2
        Transform((y.Pow(2), x.Pow(3) - 7 * x - 6), TorsionMeth.NagellLutz);
        Transform((y.Pow(2), x.Pow(3) - 2 * x + 1), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) + x.Pow(2)), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 6 * x - 4), TorsionMeth.NagellLutz);

        // Torsion C5
        Transform((y.Pow(2) + y, x.Pow(3) - x.Pow(2)), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) + 15 * x + 9), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) + x.Pow(2) + 1), TorsionMeth.NagellLutz);

        // Torsion C6
        Transform((y.Pow(2), x.Pow(3) + 1), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - x.Pow(2) + 6 * x), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - 6 * x + 4), TorsionMeth.NagellLutz);

        // Torsion C7
        Transform((y.Pow(2) + x * y, x.Pow(3) + 159 * x + 1737), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - x + 137), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 3 * x + 3), TorsionMeth.NagellLutz);

        // Torsion C8, C4 x C2
        Transform((y.Pow(2), x.Pow(3) + x.Pow(2) + 16 * x + 180), TorsionMeth.NagellLutz);
        Transform((y.Pow(2), x.Pow(3) - x.Pow(2) - 4 * x + 4), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - 34 * x + 68), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - 4 * x - 1), TorsionMeth.NagellLutz);

        // Torsion C9
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 14 * x + 29), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) + 108 * x + 11664), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - 4767 * x + 127449), TorsionMeth.NagellLutz);

        // Torsion C10
        Transform((y.Pow(2) + x * y, x.Pow(3) - 45 * x + 81), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) + 115 * x + 561), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - 828 * x + 9072), TorsionMeth.NagellLutz);

        // Torsion C12, C6 x C2
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - 19 * x + 26), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 122 * x + 1721), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - 361 * x + 2585), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y + y, x.Pow(3) + 1922 * x + 20756), TorsionMeth.NagellLutz);

        // Torsion C8 x C2
        Transform((y.Pow(2) + x * y, x.Pow(3) - 1070 * x + 7812), TorsionMeth.NagellLutz);
        Transform((y.Pow(2) + x * y, x.Pow(3) - 8696090 * x + "9838496100"), TorsionMeth.NagellLutz);

        GlobalStopWatch.Show("END"); // Time:5.955s
    }

    public static void Example6MordellCurve()
    {
        var listMissing = new List<int>();
        for (int n = 1; n <= 100; n++)
        {
            var (gEll, ellpts, abType) = NagellLutzTorsionGroup([0, n], LogLevel.Off);
            var E = (EllGroup<Rational>)gEll.BaseGroup;
            var intPts = ellpts.ToHashSet();

            var sz = -1;
            while (sz != intPts.Count)
            {
                sz = intPts.Count;
                foreach (var pt in intPts.ToArray())
                    AddNewIntegralPoint(E, intPts, pt);
            }

            var missing = GroupExt.B081119[n] - intPts.Count(pt => !pt.IsO);
            if (missing != 0)
            {
                listMissing.Add(n);
                Console.WriteLine($"## y^2 = x^3 + {n}");
                if (abType.Length == 1 && abType[0] == 1)
                    Console.WriteLine("Trivial Torsion. TODO Integer Points"); // TODO: Trivial torsion

                Console.WriteLine($"Sols {{ {intPts.Glue(", ")} }} A081119 missing {missing}");
            }
            else if (intPts.Count(pt => !pt.IsO) != 0)
                intPts.Where(pt => !pt.IsO).Println($"{intPts.Count(pt => !pt.IsO)} Integer Points of y^2 = x^3 + {n}");

            Console.WriteLine();
        }

        Console.WriteLine($"listMissing = {listMissing.Count} [{listMissing.Glue(", ")}]");
    }

    public static void Example7MordellsEquation()
    {
        GlobalStopWatch.Restart();
        var mordell = new Dictionary<int, (BigInteger x, BigInteger y)[]>();
        int[][] xmax10e4 =
            [[10000], [17, 24, 100, 141, 217, 388, 414, 513, 516, 521, 568, 649, 659, 740, 757, 836, 960, 985]];
        int[][] xmax10e5 = [[100000], [297, 377, 427, 885, 899]];
        int[][] xmax10e6 = [[1000000], [225, 353, 618]];
        var list = new List<int[][]>() { xmax10e4, xmax10e5, xmax10e6 };
        var nmax = 1000;
        var tors = new Dictionary<int, int[]>();
        foreach (int n in nmax.Range(1))
        {
            var xmin = (int)double.Ceiling(double.Pow(n, 1.0 / 3.0));
            var idx = list.FindIndex(e => e[1].Contains(n));
            var xmax = idx != -1 ? list[idx][0][0] : 1000;
            var set = new HashSet<(BigInteger, BigInteger)>();
            for (var x = -xmin - 1; x < xmax; x++)
            {
                var (y, info) = ApproxSolver(x, n);
                if (!info)
                    continue;

                set.UnionWith([(x, y), (x, -y)]);
            }

            var pts = mordell[n] = set.Order().ToArray();
            var abType = NagellLutzTorsionGroup([0, n], pts, LogLevel.Off).abType;
            // var abType = NagellLutz(0, n, pts.Select(e => new EllPt<Rational>($"{e.x}", $"{e.y}")).ToArray()).abType;
            tors[n] = abType;
            Console.WriteLine($"n = {n} ");
            Console.CursorTop--;
        }

        var missingSet = new List<int>();
        foreach (var (n, pts) in mordell)
        {
            var nb = pts.Length;
            var missing = GroupExt.B081119[n] - nb;
            Console.WriteLine($"n = {n,4}, Torsion {tors[n].Glue(" x ", "C{0}")}," +
                              $" Integral Points {nb,-3}/{GroupExt.B081119[n],3} => {{ {pts.Glue(", ")} }}");
            if (pts.Any(e => BigInteger.Pow(e.y, 2) != BigInteger.Pow(e.x, 3) + n))
                throw new();

            if (missing != 0)
                missingSet.Add(n);
        }

        if (missingSet.Count != 0)
            Console.WriteLine($"{missingSet.Count} Missing {{ {missingSet.Glue(", ")} }}");

        Console.WriteLine();
        GlobalStopWatch.Show(); // Time:5.009s
    }

    // Daniel Guin - Thomas Hausberger, Algebre Tome 1, page 174
    public static void Example8ConguentNumbers()
    {
        var primXYZ = 10.Range(1).Grid2D().Select(e => (a: e.t1, b: e.t2))
            .Where(e => e.a < e.b && (e.a * e.b) % 2 == 0 && IntExt.Gcd(e.a, e.b) == 1)
            .Select(e => (x0: e.b * e.b - e.a * e.a, y0: 2 * e.a * e.b, z: e.a * e.a + e.b * e.b))
            .Select(e => (x: BigInteger.Min(e.x0, e.y0), y: BigInteger.Max(e.x0, e.y0), z: new BigInteger(e.z)))
            .Distinct().Order().ToArray();

        primXYZ.Println(e => $"{e.z}^2 = {e.x}^2 + {e.y}^2", $"Z^2 = X^2 + Y^2, Count:{primXYZ.Length}");
        Console.WriteLine();

        var lt = new HashSet<(Rational X, Rational Y, Rational Z)>();
        foreach (var (x, y, z) in primXYZ)
        {
            var n = (x * y) / 2;
            Console.WriteLine($"{z}^2 = {x}^2 + {y}^2 and n = {n}");

            var n2 = new Rational(n).Pow(2);
            EllRank(-n * n, 0);
            var E = new EllGroup<Rational>(-n2, "0");
            var (gEll, intPts, _) = NagellLutzTorsionGroup(E, LogLevel.Off);
            foreach (var pt in intPts.ToArray())
                AddNewIntegralPoint(E, intPts, pt);

            Console.WriteLine($"n = {n} Ell y^2 = x^3 - {n2}x");
            var pts = intPts.Where(e => e.IsIntegral() && !e.IsO && !e.Y.IsZero())
                .Select(e => (pt: e, X: (n2 - e.X.Pow(2)) / e.Y, Y: -2 * n * e.X / e.Y, Z: (n2 + e.X.Pow(2)) / e.Y))
                .Select(e => (e.pt, X: e.X.Absolute, Y: e.Y.Absolute, Z: e.Z.Absolute))
                .DistinctBy(e => (Rational.Min(e.X, e.Y), Rational.Max(e.X, e.Y)))
                .ToArray();

            lt.UnionWith(pts.Select(e => (Rational.Min(e.X, e.Y), Rational.Max(e.X, e.Y), e.Z)));

            pts.Println(e => $"{e.pt,-20} X:{e.X} Y:{e.Y} Z:{e.Z}", $"Pts = {pts.Length}");
            Console.WriteLine($"Check points X^2+Y^2=Z^2:{pts.All(e => (e.X.Pow(2) + e.Y.Pow(2)).Equals(e.Z.Pow(2)))}");
            Console.WriteLine();
        }

        lt.OrderBy(e => e.Z).ThenBy(e => e.X)
            .Println(e => $"X:{e.X} Y:{e.Y} Z:{e.Z}", $"All X^2+Y^2=Z^2 total {lt.Count}");
    }

    public static void Example9Rank()
    {
        GlobalStopWatch.Restart();

        // Rank 0
        EllRank(-432, 8208);
        EllRank(-675, 13662);
        EllRank(-27, 8694);

        // Rank 1
        EllRank(0, 3);
        EllRank(-36, 0);
        EllRank(-961, 0);

        // Rank 2
        EllRank(-3024, 46224);
        EllRank(-5292, -101520);
        EllRank(-2052, 34560);

        // Rank 0, 1, 2, 3, 4
        foreach (var d in new[] { 1, 5, 34, 1254, 29274 })
            EllRank(-d * d, 0);

        GlobalStopWatch.Show(); // Time:5.079s
    }
    
    public static void Example10ECDLP()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        Logger.Level = LogLevel.Level1;
    
        // C : y^2 = x^3 + x^2 + x + 1 p = 97, P = (7, 20), Q = (17, 46)
        // sage: E=EllipticCurve(GF(97),[0,1,0,1,1]);P=E([7, 20]);Q=Ea([17, 46]);Q.log(P)
        // pari: E=ellinit([0,1,0,1,1],97);elllog(E,[17, 46],[7, 20])
        ECDLP(p:97, (7, 20), (17, 46), [0, 1, 0, 1, 1]);
    
        // (a) C : y^2 = x^3 + x^2 + x + 3, p = 103, P = (7, 14), Q = (8, 22).
        ECDLP(p:103, (7, 14), (8, 22), [0, 1, 0, 1, 3]);

        // (b) C : y^2 = x^3 − 2x^2 + 5x + 6, p = 149, P = (11, 16), Q = (110, 46).
        ECDLP(p:149, (11, 16), (110, 46), [0, -2, 0, 5, 6]);

        // (c) C : y^2 = x^3 + x^2 + x + 2, p = 10037, P = (8, 7358), Q = (2057, 5437).
        ECDLP(p:10037, (8, 7358), (2057, 5437), [0, 1, 0, 1, 2]);

        // C : y^2 + x*y + y = x^3 - x^2 - 29*x - 53, p = 10037, P = (8, 7358), Q = (2057, 5437).
        ECDLP(p: 20011, (16897, 9208), (8965, 18468), [1, -1, 1, -29, -53]);
        // |Ell[1,20010,1,19982,19958](Z/20011Z)| = 19872
        // P=(16897, 9208) Q=( 8965,18468) 4305xP=Q
        // #  Time:3ms
    }
}