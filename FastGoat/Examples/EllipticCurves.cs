using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class EllipticCurves
{
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

    static BigCplx[] NRoots(KPoly<BigCplx> P, int maxLoop = 200)
    {
        var O1 = P.KZero.O;
        var O2 = O1 / 2;
        var P0 = P;
        var roots = new List<BigCplx>();
        var (pi, e) = (BigReal.Pi(O2), BigReal.E(O2));
        var a0 = new BigCplx(pi, e);

        while (P0.Degree > 0)
        {
            a0 = new BigCplx(pi, e);
            var a1 = FG.NSolve(P0, a0.ToBigCplx(O1), maxLoop);
            roots.Add(a1);
            P0 /= P0.X - a1;
        }

        return roots.ToArray();
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
        Console.WriteLine($"Ell[{a},{b}](Q) Torsion = {tor!.Descending().Glue(" x ", "C{0}")}");
        Console.WriteLine();
    }

    static (BigInteger disc, BigInteger[]) CandidatsY(BigInteger a, BigInteger b)
    {
        var disc = 4 * BigInteger.Pow(a, 3) + 27 * BigInteger.Pow(b, 2);
        Console.WriteLine(new { disc });
        var r = IntExt.PrimesDec(BigInteger.Abs(disc)).Where(e => e.Value >= 2)
            .Aggregate(BigInteger.One, (acc, e) => acc * BigInteger.Pow(e.Key, e.Value / 2));

        var divs = IntExt.DividorsBigInt(r).Order().ToArray();
        return (disc, divs);
    }

    static IEnumerable<EllPt<Rational>> SolveX(BigInteger a, BigInteger b, BigInteger[] Ys, bool show = false)
    {
        var x = FG.BCplxPoly();
        var O = x.KZero.O;
        var (A, B) = (BigCplx.FromRational(new(a), O), BigCplx.FromRational(new(b), O));
        foreach (var y0 in Ys.Prepend(0))
        {
            var y = BigCplx.FromRational(new(y0), O);
            var P = x.Pow(3) + A * x + B - y.Pow(2);
            var sols = NRoots(P);
            if (show)
                sols.Println($"Y = {y} P = {P}");

            var ellpts = sols.Where(xi => xi.IsInteger())
                .Select(xi => new EllPt<Rational>(xi.RealPart.Round0.ToRational, y.RealPart.ToRational));

            foreach (var pt in ellpts)
                yield return pt;
        }
    }

    static void NagellLutz(BigInteger a, BigInteger b, bool show = false)
    {
        var (disc, Ys) = CandidatsY(a, b);
        var ellpts = SolveX(a, b, Ys, show).ToArray();
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
            }

            if (acc.IsO)
                set.Add(pt);
        }

        if (show)
            set.Println("Elements");

        var gEll = Group.Generate(E, set.ToArray());
        DisplayGroup.HeadElements(gEll);
        var abType = Group.AbelianGroupType(gEll);
        Console.WriteLine($"{gEll} Torsion = {abType.Glue(" x ", "C{0}")}");
        Console.WriteLine();
    }

    public static void Transform((Polynomial<Rational, Xi> lhs, Polynomial<Rational, Xi> rhs) e, bool nagellLutz = true)
    {
        Console.WriteLine($"{e.lhs} = {e.rhs}");
        var F = -e.lhs + e.rhs;
        var ((y, Y), (x, X)) = F.IndeterminatesAndVariables.Deconstruct();

        var a = Ring.Decompose(F, y).Item1[Y].ConstTerm;
        Console.WriteLine($"{F} = 0");
        if (!a.IsZero())
        {
            Console.WriteLine($"{Y} <- {Y + a / 2}");
            F = F.Substitute(Y + a / 2, y);
            Console.WriteLine($"{F} = 0");
        }

        var cX = Ring.Decompose(Ring.Decompose(F, x).Item1[X], y).Item1;
        var b = cX.ContainsKey(Y) ? cX[Y] : F.Zero;
        if (!b.IsZero())
        {
            Console.WriteLine($"{Y} <- {Y + b * X / 2}");
            F = F.Substitute(Y + b * X / 2, y);
            Console.WriteLine($"{F} = 0");
        }

        var c = Ring.Decompose(F, x).Item1[X.Pow(2)];
        if (!c.IsZero())
        {
            Console.WriteLine($"{X} <- {X - c / 3}");
            F = F.Substitute(X - c / 3, x);
            Console.WriteLine($"{F} = 0");
        }

        Console.WriteLine($"{Y.Pow(2)} = {F + Y.Pow(2)}");
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
            F = (d1.Pow(3) * F.Substitute(X / d1, x)).Substitute(Y / d2, y);
            Console.WriteLine($"{F} = 0");
            Console.WriteLine();
        }

        var P1 = (F + Y.Pow(2)).ToKPoly(x);
        var (A, B) = (P1[1], P1[0]);
        EllTors(A.Num, B.Num);
        if (nagellLutz)
            NagellLutz(A.Num, B.Num);
        Console.WriteLine();
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
        EllTors(-43, 166, nbPrimes: 20); // Ell[-43,166](Q) Torsion = C7
    }

    public static void ExampleNagellLutz()
    {
        NagellLutz(-36, 0, show: true); // Ell[-36,0](Q) Torsion = C2 x C2
        NagellLutz(0, 3); // Ell[0,3](Q) Torsion = C1
        NagellLutz(1, 0); // Ell[1,0](Q) Torsion = C2
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

        Transform(e1);
        Transform(e2);
        Transform(e3);
        Transform(e4);
        Transform(e5);
        Transform(e6);
        Transform(e7);
        Transform(e8);
    }
}