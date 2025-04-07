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

    static int[] EllFp(BigInteger a, BigInteger b, int p, bool show = false)
    {
        var (A, B) = (new ZnBigInt(p, a), new ZnBigInt(p, b));
        var disc = 4 * A.Pow(3) + 27 * B.Pow(2);
        if (disc.IsZero())
            return [];

        var E = new EllGroup<ZnBigInt>(A, B);
        var ell = p.Range().Select(k => new ZnBigInt(p, k))
            .Select(x => (x, y2: x.Pow(3) + A * x + B))
            .Select(e => (e.x, y: NumberTheory.SqrtModANTV1(e.y2.K, p) * e.x.One))
            .Select(e => new EllPt<ZnBigInt>(e.x, e.y))
            .Where(e => E.Contains(e.X, e.Y))
            .Order()
            .ToArray();

        var gEll = Group.Generate(E, ell);
        if (show)
            DisplayGroup.HeadElements(gEll);

        var abType = Group.AbelianGroupType(gEll);
        Console.WriteLine($"{gEll} ~ {abType.Glue(" x ", "C{0}")}");
        if (gEll.Any(e => !e.IsO && !E.Contains(e.X, e.Y)))
            throw new();

        return abType;
    }

    static void EllTors(BigInteger a, BigInteger b, int nbPrimes = 10, bool show = false)
    {
        Console.WriteLine($"#### Start Ell[{a},{b}](Q)");
        var disc = 4 * BigInteger.Pow(a, 3) + 27 * BigInteger.Pow(b, 2);
        var allTypes = IntExt.Primes10000.Where(p => (2 * disc) % p != 0).Take(nbPrimes)
            .Select(n => EllFp(a, b, n, show))
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

        allTypes.Println(e => e.Glue(" x ", "C{0}"), $"Morphisms Ell[{a},{b}](Q) ->");
        set.Select(l => l.ToArray()).Println(e => e.Glue(" x ", "C{0}"), "Intersections subgroups");
        var tor = set.MaxBy(l => l.Aggregate(1, (acc, i) => acc * i));
        Console.WriteLine($"Ell[{a},{b}](Q) Torsion ~-> {tor!.Descending().Glue(" x ", "C{0}")}");
        Console.WriteLine();
    }

    static (BigInteger disc, BigInteger[]) CandidatsY(BigInteger a, BigInteger b)
    {
        var disc = 16 * (4 * BigInteger.Pow(a, 3) + 27 * BigInteger.Pow(b, 2));
        Console.WriteLine(new { disc });
        var r = IntExt.PrimesDec(BigInteger.Abs(disc))
            .Aggregate(BigInteger.One, (acc, r) => acc * BigInteger.Pow(r.Key, r.Value / 2 + r.Value % 2));
        var divs = IntExt.DividorsBigInt(r).Where(y => disc % (y * y) == 0).Order().ToArray();
        return (disc, divs);
    }

    static IEnumerable<EllPt<Rational>> SolveX(BigInteger a, BigInteger b, BigInteger[] Ys, bool show = false)
    {
        var X = FG.QPoly();
        var (A, B) = (new Rational(a), new Rational(b));
        foreach (var y in Ys.Prepend(0).Select(y => new Rational(y)))
        {
            var P = X.Pow(3) + A * X + B;
            var sols = IntFactorisation.FactorsQ(P - y.Pow(2));
            if (show)
                sols.Println($"Y = {y}, solve {y.Pow(2)} = {P}");

            var ellpts = sols.Where(e => e.Item1.Degree == 1).Select(e => new EllPt<Rational>(-e.Item1[0], y));
            foreach (var pt in ellpts)
                yield return pt;
        }
    }

    static (EllGroup<Rational> E, ConcreteGroup<EllPt<Rational>> gEll, int[] abType, HashSet<EllPt<Rational>> intPts)
        NagellLutz(BigInteger a, BigInteger b, EllPt<Rational>[] ellpts)
    {
        var intPts = ellpts.Concat(ellpts.Select(e => new EllPt<Rational>(e.X, -e.Y))).ToHashSet();
        var E = new EllGroup<Rational>($"{a}", $"{b}");
        var set = new List<EllPt<Rational>>() { E.O };
        foreach (var pt in ellpts)
        {
            var acc = E.O;
            for (int i = 1; i <= 12; i++)
            {
                acc = E.Op(acc, pt);
                if (acc.IsO || !acc.X.IsInteger() || !acc.Y.IsInteger())
                    break;

                intPts.Add(acc);
                intPts.Add(new(acc.X, -acc.Y));
            }

            if (acc.IsO)
                set.Add(pt);
        }

        var gEll = Group.Generate(E, set.ToArray());
        var abType = Group.AbelianGroupType(gEll);
        return (E, gEll, abType, intPts);
    }

    static (int[] abType, HashSet<EllPt<Rational>> intPts, EllGroup<Rational> E)
        NagellLutz(BigInteger a, BigInteger b, Func<EllPt<Rational>, EllPt<Rational>> revTrans, bool show = false)
    {
        var (disc, Ys) = CandidatsY(a, b);
        var ellpts = SolveX(a, b, Ys, show).ToArray();
        var (E, gEll, abType, intPts) = NagellLutz(a, b, ellpts);

        if (show)
            intPts.Println("Elements");

        DisplayGroup.HeadElements(gEll);
        Console.WriteLine($"{gEll} Torsion = {abType.Glue(" x ", "C{0}")}");
        Console.WriteLine();
        var intPtsF = gEll.Concat(intPts).Select(pt => revTrans(pt)).ToHashSet();
        return (abType, intPtsF, E);
    }

    static (int[] abType, HashSet<EllPt<Rational>> intPts, EllGroup<Rational> E)
        NagellLutz(BigInteger a, BigInteger b, bool show = false)
    {
        return NagellLutz(a, b, pt => pt, show);
    }

    public static (EllGroup<Rational> E, ConcreteGroup<EllPt<Rational>> gEll, int[] abType,
        HashSet<EllPt<Rational>> intPts) NagellLutzTorsionGroup(BigInteger a, BigInteger b)
    {
        var (disc, Ys) = CandidatsY(a, b);
        var ellpts = SolveX(a, b, Ys, show: false).ToArray();
        return NagellLutz(a, b, ellpts);
    }

    static (KPoly<Rational> P1, Func<EllPt<Rational>, EllPt<Rational>> revTrans)
        MinimizedForm((Polynomial<Rational, Xi> lhs, Polynomial<Rational, Xi> rhs) e)
    {
        var F = -e.lhs + e.rhs;
        var ((y, Y), (x, X)) = F.IndeterminatesAndVariables.Deconstruct();
        var ind = X.Indeterminates;
        var (xm, ym) = (new Monom<Xi>(ind, x), new Monom<Xi>(ind, y));
        var (xym, x2m) = (xm.Mul(ym), xm.Pow(2));
        var (a1, a2, a3, a4, a5) = (F[xym], F[ym], F[x2m], F[xm], F.ConstTerm);
        var A = -a1.Pow(4) / 48 - a1.Pow(2) * a3 / 6 + a1 * a2 / 2 - a3.Pow(2) / 3 + a4;
        var B = a1.Pow(6) / 864 + a1.Pow(4) * a3 / 72 - a1.Pow(3) * a2 / 24 + a1.Pow(2) * a3.Pow(2) / 18 -
            a1.Pow(2) * a4 / 12 - a1 * a2 * a3 / 6 + 2 * a3.Pow(3) / 27 + a2.Pow(2) / 4 - a3 * a4 / 3 + a5;

        var sqDivs864 = new[] { 1, 4, 9, 16, 36, 144 } // Square Divisidors of 864 
            .Select(div => (div, pow2: div * div, pow3: div * div * div)).ToArray();
        var (sqDiv, _, sqDivPow3) = sqDivs864.OrderBy(f => f.div)
            .First(div => div.pow2 % A.Denom == 0 && div.pow3 % B.Denom == 0);
        var d1 = new Rational(sqDiv);
        var d2 = new Rational(IntExt.SqrtBigInt(sqDivPow3));

        Func<EllPt<Rational>, EllPt<Rational>> revTrans = pt =>
        {
            if (pt.IsO)
                return pt;
        
            var _x = pt.X / d1 - a1.Pow(2) / 12 - a3 / 3;
            var _y = pt.Y / d2;
            return new(_x, -_y + (a1 * _x + a2) / 2);
        };

        A *= d1.Pow(2);
        B *= d1.Pow(3);
        if (Logger.Level != LogLevel.Off)
        {
            Console.WriteLine($"Elliptic curve      {e.lhs} = {e.rhs}");
            Console.WriteLine($"Simplified form     y^2 = {X.Pow(3) + A * X + B}");
            Console.WriteLine();
        }

        var P1 = (X.Pow(3) + A * X + B).ToKPoly(x);
        return (P1, revTrans);
    }

    static void Transform((Polynomial<Rational, Xi> lhs, Polynomial<Rational, Xi> rhs) e,
        TorsionMeth meth = TorsionMeth.Both)
    {
        var (P1, revTrans) = MinimizedForm(e);
        Console.WriteLine($"Elliptic curve      {e.lhs} = {e.rhs}");
        Console.WriteLine($"Simplified form     y^2 = {P1}");
        var (A, B) = (P1[1], P1[0]);

        if ((meth & TorsionMeth.Fp) == TorsionMeth.Fp)
            EllTors(A.Num, B.Num);
        if ((meth & TorsionMeth.NagellLutz) == TorsionMeth.NagellLutz)
        {
            var ng = NagellLutz(A.Num, B.Num, revTrans);
            ng.intPts.OrderBy(pt => pt.Height())
                .Println($"Points of Elliptic curve      {e.lhs} = {e.rhs}");
            Console.WriteLine();
        }
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

    static int SchoofEllPtsCount(BigInteger a, BigInteger b, int p)
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
        EllTors(-36, 0, show: true); // Ell[-36,0](Q) Torsion = C2 x C2
        EllTors(0, 3); // Ell[0,3](Q) Torsion = C1
        EllTors(1, 0); // Ell[1,0](Q) Torsion = C2
        EllTors(0, 1); // Ell[0,1](Q) Torsion = C6
        EllTors(-43, 166, nbPrimes: 20); // Ell[-43,166](Q) Torsion = C7
    }

    public static void Example3NagellLutz()
    {
        NagellLutz(-36, 0, show: true); // Ell[-36,0](Q) Torsion = C2 x C2
        NagellLutz(0, 3); // Ell[0,3](Q) Torsion = C1
        NagellLutz(1, 0); // Ell[1,0](Q) Torsion = C2
        NagellLutz(0, 1); // Ell[0,1](Q) Torsion = C6
        NagellLutz(-43, 166); // Ell[-43,166](Q) Torsion = C7
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

        GlobalStopWatch.Show("END"); // Time:35.122s
    }

    public static void Example6MordellCurve()
    {
        for (int n = 1; n <= 100; n++)
        {
            var (abType, ellpts, E) = NagellLutz(0, n);
            var intPts = ellpts.ToHashSet();
            var sz = -1;
            while (sz != intPts.Count)
            {
                sz = intPts.Count;
                intPts.UnionWith(intPts.Prepend(E.O).ToArray().Grid2D()
                    .SelectMany(e => new[] { E.Op(e.t1, e.t2), E.Invert(E.Op(e.t1, e.t2)) })
                    .Where(pt => !pt.IsO && pt.X.IsInteger() && pt.Y.IsInteger()).ToArray());
            }

            var missing = GroupExt.B081119[n] - intPts.Count(pt => !pt.IsO);
            if (missing != 0)
            {
                Console.WriteLine($"## y^2 = x^3 + {n}");
                if (abType.Length == 1 && abType[0] == 1)
                    Console.WriteLine("Trivial Torsion. TODO Integer Points"); // TODO: Trivial torsion

                Console.WriteLine($"Sols {{ {intPts.Glue(", ")} }} A081119 missing {missing}");
            }
            else if (intPts.Count(pt => !pt.IsO) != 0)
                intPts.Where(pt => !pt.IsO).Println($"{intPts.Count(pt => !pt.IsO)} Integer Points of y^2 = x^3 + {n}");

            Console.WriteLine();
        }
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
            var abType = NagellLutz(0, n, pts.Select(e => new EllPt<Rational>($"{e.x}", $"{e.y}")).ToArray()).abType;
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
            .Select(e => (x: int.Min(e.x0, e.y0), y: int.Max(e.x0, e.y0), e.z))
            .Order()
            .ToArray();

        primXYZ.Println(e => $"{e.z}^2 = {e.x}^2 + {e.y}^2", $"Z^2 = X^2 + Y^2, Count:{primXYZ.Length}");
        Console.WriteLine();
        // primXYZ.Where(e => e.x < 65000 && e.y < 65000 && e.x * e.y > 0).Select(e => e.x * e.y / 2).Select(e =>
        //         IntExt.PrimesDec(e).Select(f => (f.Key, f.Value % 2))
        //             .Aggregate(1, (acc, f) => acc * f.Key.Pow(f.Item2)))
        //     .Distinct().Order().Take(1000).Println($"n =");
        // Console.ReadLine();

        foreach (var (x, y, z) in primXYZ)
        {
            var n = (x * y) / 2;
            Console.WriteLine($"{z}^2 = {x}^2 + {y}^2 and n = {x * y / 2}");

            var (a, b) = (-n * n, 0);
            var (disc, Ys) = CandidatsY(a, b);
            var ellpts = SolveX(a, b, Ys).ToArray();
            EllRank(a, b);
            var (E, gEll, abType, intPts) = NagellLutz(a, b, ellpts);

            Console.WriteLine($"n = {n} Ell y^2 = x^3 - {n * n}x");
            var pts = intPts.Where(e => !e.IsO && !e.Y.IsZero()).Select(e => (pt: e,
                X: (n.Pow(2) - e.X.Pow(2)) / e.Y, Y: -2 * n * e.X / e.Y, Z: (n.Pow(2) + e.X.Pow(2)) / e.Y)).ToArray();

            pts.Println(e => $"{e.pt,-20} X:{e.X} Y:{e.Y} Z:{e.Z}", "Pts");
            Console.WriteLine($"Check points X^2+Y^2=Z^2:{pts.All(e => (e.X.Pow(2) + e.Y.Pow(2)).Equals(e.Z.Pow(2)))}");
            Console.WriteLine();
        }
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
}