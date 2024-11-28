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

    static int[][] AbSubTypes(int[] type)
    {
        var all = type.Select(t => IntExt.Dividors(t).Append(t).ToArray()).MultiLoop()
            .Select(l => l.Order())
            .ToHashSet(new SequenceEquality<int>())
            .Select(l => l.ToArray())
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
        var (A, B) = (new ZnInt(p, (int)BigInteger.Remainder(a, p)), new ZnInt(p, (int)BigInteger.Remainder(b, p)));
        var disc = 4 * A.Pow(3) + 27 * B.Pow(2);
        if (disc.IsZero())
            return [];

        var E = new EllGroup<ZnInt>(A, B);
        var ell = p.Range().Select(k => new ZnInt(p, k)).ToArray().Grid2D()
            .Select(e => new EllPt<ZnInt>(e.t1, e.t2))
            .Where(e => E.Contains(e.X, e.Y))
            .Distinct()
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

        allTypes.Println(e => e.Glue(" x ", "C{0}"), $"Morphism Ell[{a},{b}](Q) ->");
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
        var divs = IntExt.DividorsBigInt(r).Order().ToArray();
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

    static (int[] abType, HashSet<EllPt<Rational>> intPtsF, EllGroup<Rational> E) NagellLutz(BigInteger a, BigInteger b,
        List<Func<EllPt<Rational>, EllPt<Rational>>> revTrans, bool show = false)
    {
        var (disc, Ys) = CandidatsY(a, b);
        var ellpts = SolveX(a, b, Ys, show).ToArray();
        var (E, gEll, abType, intPts) = NagellLutz(a, b, ellpts);

        if (show)
            intPts.Println("Elements");

        DisplayGroup.HeadElements(gEll);
        Console.WriteLine($"{gEll} Torsion = {abType.Glue(" x ", "C{0}")}");
        Console.WriteLine();
        revTrans.Reverse();
        var intPtsF = gEll.Where(pt => !pt.IsO).Select(pt => revTrans.Aggregate(pt, (acc, trans) => trans(acc)))
            .Concat(intPts)
            .ToHashSet();

        return (abType, intPtsF, E);
    }

    static (int[] abType, HashSet<EllPt<Rational>> intPtsF, EllGroup<Rational> E)
        NagellLutz(BigInteger a, BigInteger b, bool show = false)
    {
        return NagellLutz(a, b, new(), show);
    }

    static void Transform((Polynomial<Rational, Xi> lhs, Polynomial<Rational, Xi> rhs) e,
        TorsionMeth meth = TorsionMeth.Both)
    {
        Console.WriteLine($"Initial   form {e.lhs} = {e.rhs}");
        var F = -e.lhs + e.rhs;
        var ((y, Y), (x, X)) = F.IndeterminatesAndVariables.Deconstruct();

        var a = Ring.Decompose(F, y).Item1[Y].ConstTerm;
        Console.WriteLine($"{F} = 0");
        var revTrans = new List<Func<EllPt<Rational>, EllPt<Rational>>>();

        if (!a.IsZero())
        {
            Console.WriteLine($"{Y} <- {Y + a / 2}");
            F = F.Substitute(Y + a / 2, y);
            revTrans.Add(pt => new(pt.X, pt.Y + a / 2));
            Console.WriteLine($"{F} = 0");
        }

        var cX = Ring.Decompose(Ring.Decompose(F, x).Item1[X], y).Item1;
        var b = cX.ContainsKey(Y) ? cX[Y] : F.Zero;
        if (!b.IsZero())
        {
            Console.WriteLine($"{Y} <- {Y + b * X / 2}");
            F = F.Substitute(Y + b * X / 2, y);
            revTrans.Add(pt => new(pt.X, pt.Y + b.ConstTerm * pt.X / 2));
            Console.WriteLine($"{F} = 0");
        }

        var c = Ring.Decompose(F, x).Item1[X.Pow(2)];
        if (!c.IsZero())
        {
            Console.WriteLine($"{X} <- {X - c / 3}");
            F = F.Substitute(X - c / 3, x);
            revTrans.Add(pt => new(pt.X - c.ConstTerm / 3, pt.Y));
            Console.WriteLine($"{F} = 0");
        }

        var P = F + Y.Pow(2);
        if (P.NbIndeterminates != 1)
            throw new();

        var P0 = P.ToKPoly(x);
        var (a0, b0) = (P0[1].Denom, P0[0].Denom);
        var decompA = !a0.IsZero
            ? IntExt.PrimesDec(a0).ToDictionary(r => r.Key, r => (r.Value % 2 == 0 ? 0 : 1) + r.Value / 2)
            : new();
        var decompB = !b0.IsZero
            ? IntExt.PrimesDec(b0).ToDictionary(r => r.Key, r => (r.Value % 3 == 0 ? 0 : 1) + r.Value / 3)
            : new();
        foreach (var i in decompA.Keys.Intersect(decompB.Keys))
        {
            var (ai, bi) = (decompA[i], decompB[i]);
            decompB[i] = int.Min(ai, bi);
            decompA.Remove(i);
        }

        var decompD = decompA.Concat(decompB).ToDictionary(r => r.Key, r => r.Value + r.Value % 2);
        var d1 = new Rational(decompD.Aggregate(BigInteger.One, (acc, r) => acc * r.Key.Pow(r.Value)));
        var d2 = new Rational(decompD.Aggregate(BigInteger.One, (acc, r) => acc * r.Key.Pow(3 * r.Value / 2)));
        if (!(d1 - 1).IsZero())
        {
            Console.WriteLine($"{X} <- {X / d1}");
            Console.WriteLine($"{Y} <- {Y / d2}");
            revTrans.Add(pt => new(pt.X / d1, pt.Y / d2));
            F = (d1.Pow(3) * F.Substitute(X / d1, x)).Substitute(Y / d2, y);
            Console.WriteLine($"{F} = 0");
        }

        Console.WriteLine($"Minimized form {Y.Pow(2)} = {F + Y.Pow(2)}");
        Console.WriteLine();

        var P1 = (F + Y.Pow(2)).ToKPoly(x);
        var (A, B) = (P1[1], P1[0]);

        if ((meth & TorsionMeth.Fp) == TorsionMeth.Fp)
            EllTors(A.Num, B.Num);
        if ((meth & TorsionMeth.NagellLutz) == TorsionMeth.NagellLutz)
        {
            var sols = NagellLutz(A.Num, B.Num, revTrans);
            // TODO: fix integral points
            // sols.Where(pt => pt.X.IsInteger() && pt.Y.IsInteger()).Order()
            //     .Println($"Integral Points of {e.lhs} = {e.rhs}"); 
        }

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
    public static void ExampleFp()
    {
        EllTors(-36, 0, show: true); // Ell[-36,0](Q) Torsion = C2 x C2
        EllTors(0, 3); // Ell[0,3](Q) Torsion = C1
        EllTors(1, 0); // Ell[1,0](Q) Torsion = C2
        EllTors(0, 1); // Ell[0,1](Q) Torsion = C6
        EllTors(-43, 166, nbPrimes: 20); // Ell[-43,166](Q) Torsion = C7
    }

    public static void ExampleNagellLutz()
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

    public static void ExampleTransformCurve()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var (x, y) = Ring.Polynomial(Rational.KZero(), "x", "y").Deconstruct();

        // – E1 : y 2 + 7xy = x3 + 16x ;
        var e1 = (y.Pow(2) + 7 * x * y, x.Pow(3) + 16 * x);
        // – E2 : y 2 + xy − 5y = x3 − 5x2 ;
        var e2 = (y.Pow(2) + x * y - 5 * y, x.Pow(3) - 5 * x.Pow(2));
        // – E3 : y 2 − y = x3 − x2 ;
        var e3 = (y.Pow(2) - y, x.Pow(3) - x.Pow(2));
        // – E4 : y 2 + xy + y = x3 − x2 − 14x + 29 ;
        var e4 = (y.Pow(2) + x * y + y, x.Pow(3) - x.Pow(2) - 14 * x + 29);
        // – E5 : y 2 + xy = x3 − 45x + 81 ;
        var e5 = (y.Pow(2) + x * y, x.Pow(3) - 45 * x + 81);
        // – E6 : y 2 + 43xy − 210y = x3 − 210x2 ;
        var e6 = (y.Pow(2) + 43 * x * y - 210 * y, x.Pow(3) - 210 * x.Pow(2));
        // – E7 : y 2 + 5xy − 6y = x3 − 3x2 ;
        var e7 = (y.Pow(2) + 5 * x * y - 6 * y, x.Pow(3) - 3 * x.Pow(2));
        // – E8 : y 2 + 17xy − 120y = x3 − 60x2 .
        var e8 = (y.Pow(2) + 17 * x * y - 120 * y, x.Pow(3) - 60 * x.Pow(2));
        // – E9 : y 2 + xy = x3 − 1070x + 7812 .
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

    /* x*y + y^2 = x^3 - 1070*x + 7812
       x^3 - x*y - y^2 - 1070*x + 7812 = 0
       y <- -1/2*x + y
       x^3 + 1/4*x^2 - y^2 - 1070*x + 7812 = 0
       x <- x - 1/12
       x^3 - y^2 - 51361/48*x + 6826609/864 = 0
       y^2 = x^3 - 51361/48*x + 6826609/864
       x <- 1/36*x
       y <- 1/216*y
       x^3 - y^2 - 1386747*x + 368636886 = 0

       #### Start Ell[-1386747,368636886](Q)
       Ell[ 1, 2](Z/11Z) ~ C8 x C2
       Ell[ 2, 7](Z/13Z) ~ C8 x C2
       Ell[11,12](Z/17Z) ~ C8 x C2
       Ell[ 6, 7](Z/19Z) ~ C8 x C2
       Ell[15,16](Z/23Z) ~ C16 x C2
       Ell[ 4,22](Z/29Z) ~ C8 x C4
       Ell[ 7,14](Z/31Z) ~ C16 x C2
       Ell[13, 3](Z/37Z) ~ C8 x C4
       Ell[37,23](Z/41Z) ~ C24 x C2
       Ell[ 3,36](Z/43Z) ~ C24 x C2
       Morphism Ell[-1386747,368636886](Q) ->
           C8 x C2
           C16 x C2
           C8 x C4
           C24 x C2
       Intersections subgroups
           C1
           C2
           C4
           C8
           C2 x C2
           C2 x C4
           C2 x C8
       Ell[-1386747,368636886](Q) Torsion = C8 x C2

       { disc = -6998115764183040000 }
       |Ell[-1386747,368636886](Q)| = 16
       Type        AbelianGroup
       BaseGroup   Ell[-1386747,368636886](Q)

       Elements
       ( 1)[1] = O
       ( 2)[2] = (-1293,0)
       ( 3)[2] = (282,0)
       ( 4)[2] = (1011,0)
       ( 5)[4] = (-285,-27216)
       ( 6)[4] = (-285,27216)
       ( 7)[4] = (2307,-97200)
       ( 8)[4] = (2307,97200)
       ( 9)[8] = (-933,-29160)
       (10)[8] = (-933,29160)
       (11)[8] = (147,-12960)
       (12)[8] = (147,12960)
       (13)[8] = (1227,-22680)
       (14)[8] = (1227,22680)
       (15)[8] = (8787,-816480)
       (16)[8] = (8787,816480)

       Ell[-1386747,368636886](Q) Torsion = C8 x C2
     */

    public static void ExampleFromLMFDB()
    {
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

        GlobalStopWatch.Show("END"); // Time:29.143s
    }

    public static void ExampleMordellCurve()
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

            var missing = GroupExt.B081119[n] - intPts.Count;
            if (missing != 0)
            {
                Console.WriteLine($"## y^2 = x^3 + {n}");
                if (abType.Length == 1 && abType[0] == 1)
                    Console.WriteLine("Trivial Torsion. TODO Integer Points"); // TODO: Trivial torsion

                Console.WriteLine($"Sols {{ {intPts.Glue(", ")} }} A081119 missing {missing}");
            }
            else if (intPts.Count != 0)
                intPts.Println($"{intPts.Count} Integer Points of y^2 = x^3 + {n}");

            Console.WriteLine();
        }
    }

    public static void ExampleMordellsEquation()
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
    public static void ExampleConguentNumbers()
    {
        var primXYZ = 10.Range(1).Grid2D().Select(e => (a: e.t1, b: e.t2))
            .Where(e => e.a < e.b && (e.a * e.b) % 2 == 0 && IntExt.Gcd(e.a, e.b) == 1)
            .Select(e => (x0: e.b * e.b - e.a * e.a, y0: 2 * e.a * e.b, z: e.a * e.a + e.b * e.b))
            .Select(e => (x: int.Min(e.x0, e.y0), y: int.Max(e.x0, e.y0), e.z))
            .Order()
            .ToArray();

        primXYZ.Println(e => $"{e.z}^2 = {e.x}^2 + {e.y}^2", $"Z^2 = X^2 + Y^2, Count:{primXYZ.Length}");
        Console.WriteLine();

        foreach (var (x, y, z) in primXYZ)
        {
            var n = (x * y) / 2;
            Console.WriteLine($"{z}^2 = {x}^2 + {y}^2 and n = {x * y / 2}");
            
            var (a, b) = (-n * n, 0);
            var (disc, Ys) = CandidatsY(a, b);
            var ellpts = SolveX(a, b, Ys).ToArray();
            var (E, gEll, abType, intPts) = NagellLutz(a, b, ellpts);
            
            Console.WriteLine($"n = {n} y^2 = x^2 - {n * n}x");
            var pts = intPts.Where(e => !e.IsO && !e.Y.IsZero()).Select(e => (pt: e,
                X: (n.Pow(2) - e.X.Pow(2)) / e.Y, Y: -2 * n * e.X / e.Y, Z: (n.Pow(2) + e.X.Pow(2)) / e.Y)).ToArray();

            pts.Println(e => $"{e.pt,-20} X:{e.X} Y:{e.Y} Z:{e.Z}", "Pts");
            Console.WriteLine($"Check points X^2+Y^2=Z^2:{pts.All(e => (e.X.Pow(2) + e.Y.Pow(2)).Equals(e.Z.Pow(2)))}");
            Console.WriteLine();
        }

        {
            Rational x = "1681/49";
            Rational y = "29520/343";
            var n = 31; // TODO: find primitve Pythagorean triple X,Y,Z with n=s^2 XY/ 2 where s in Rational
            var (X, Y, Z) = ((n.Pow(2) - x.Pow(2)) / y, -2 * n * x / y, (n.Pow(2) + x.Pow(2)) / y);
            var (a, b) = (-31 * 31, 0);
            var (disc, Ys) = CandidatsY(a, b);
            var ellpts = SolveX(a, b, Ys).ToArray();
            var (E31, gEll, abType, intPts) = NagellLutz(a, b, ellpts);
            intPts.Println("Pts");
            Console.WriteLine($"({x}, {y}) is in E31 {E31.Contains(x, y)}");
            Console.WriteLine($"X:{X} Y:{Y} Z:{Z}");
            Console.WriteLine($"({Z})^2 = ({X})^2 + ({Y})^2 {(X * X + Y * Y).Equals(Z * Z)}");
        }
    }
}