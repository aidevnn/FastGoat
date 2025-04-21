using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RecomputeAllPrimesUpTo(200000);

int Ord(int p, BigInteger A) => FactorMultiplicity(p, A).mul;

bool QuadRoot(Rational a0, Rational b0, Rational c0, int p)
{
    var (a, b, c) = (a0.ToZnInt(p).K, b0.ToZnInt(p).K, c0.ToZnInt(p).K);
    if (a == 0)
        return b != 0 || c == 0;

    if (p == 2)
        return !(a == 1 && b == 1 && c == 1);

    var disc = AmodP(b * b - 4 * a * c, p);
    return disc == 0 || LegendreJacobi(disc, p) == 1;
}

int NCubicRoots(Rational b0, Rational c0, Rational d0, int p)
{
    var (b, c, d) = (b0.ToZnInt(p), c0.ToZnInt(p), d0.ToZnInt(p));
    var x = FG.ZPoly(p);
    var Q = x.Pow(3) + b * x * x + c * x + d;
    var a0 = NumberTheory.PrimitiveRootMod(p);
    var facts = IntFactorisation.FirrFsep(Q, x.KOne * a0);
    return facts.Count(f => f.g.Degree == 1);
}

TateAlgo TateAlgorithm(EllCoefs<Rational> E, int p)
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
            {
                Console.WriteLine("restart");
                return TateAlgorithm(Etmp.Transform(0, 0, 0, p), p);
            }
        }
    }
}

int EllAp(EllGroup<Rational> E, int p)
{
    if (p <= 3)
    {
        var (_a1, _a2, _a3, _a4, _a6) = E.Coefs;
        var (a1, a2, a3, a4, a6) = (_a1.ToZnInt(p), _a2.ToZnInt(p), _a3.ToZnInt(p), _a4.ToZnInt(p), _a6.ToZnInt(p));
        var card = 1 + p.Range().Select(k => new ZnInt(p, k)).ToArray().Grid2D().Select(e => (x: e.t1, y: e.t2))
            .Count(e => (e.y.Pow(2) + a1 * e.x * e.y + a3 * e.y).Equals(e.x.Pow(3) + a2 * e.x * e.x + a4 * e.x + a6));

        return p + 1 - card;
    }

    var (A, B, _, _) = E.ShortForm;
    var (a, b) = (A.ToZnBInt(new(p, 1)), B.ToZnBInt(new(p, 1)));
    return -p.Range().Select(k => new ZnBInt(p, k))
        .Select(x => (int)LegendreJacobiBigint((x.Pow(3) + a * x + b).K, p))
        .Sum(k => k <= 1 ? k : -1);
}

Dictionary<int, int> EllAn(EllGroup<Rational> E, Rational N, int maxn)
{
    var ellAn = Primes10000.Where(p => p <= maxn).ToDictionary(p => p, p => EllAp(E, p));
    ellAn[0] = 0;
    ellAn[1] = 1;
    foreach (var p in ellAn.Keys.Where(p => p > 1).Order().ToArray())
    {
        for (int i = 2; p.Pow(i) <= maxn; i++)
        {
            var t = N.Mod(p).IsZero() ? 0 : 1;
            var e1 = ellAn[p];
            var e2 = ellAn[p.Pow(i - 1)];
            var e3 = ellAn[p.Pow(i - 2)];
            var e4 = e2 * e1 - t * p * e3;
            ellAn[p.Pow(i)] = e4;
        }
    }

    for (int i = 0; i <= maxn; i++)
    {
        if (ellAn.ContainsKey(i))
            continue;

        var dec = PrimesDec(i);
        ellAn[i] = dec.Select(e => ellAn[e.Key.Pow(e.Value)]).Aggregate(1, (acc, aj) => acc * aj);
    }

    return ellAn;
}

(Rational N, Dictionary<int, int> ellAn, TateAlgo[] tate, EllGroup<Rational> Ell) EllInfos(BigInteger[] curve, int maxn)
{
    var (a1, a2, a3, a4, a6) = curve.Select(i => new Rational(i)).Deconstruct();
    var E = new EllCoefs<Rational>(a1, a2, a3, a4, a6);
    var Ell = new EllGroup<Rational>(a1, a2, a3, a4, a6);
    var dec = PrimesDec(E.Disc.Absolute.Num);
    var tate = dec.Keys.Select(p => TateAlgorithm(E, p)).ToArray();
    var N = tate.Select(e => new Rational(e.p).Pow(e.fp)).Aggregate((pi, pj) => pi * pj);
    var ellAn = EllAn(Ell, N, maxn);
    var t = double.Exp(-2 * double.Pi / double.Sqrt(N));
    var LE1 = 2 * ellAn.Where(e => e.Key != 0).Select(e => e.Value * double.Pow(t, e.Key) / e.Key).Sum();

    E.Show();
    Console.WriteLine($"Model {E.ModelStr}");
    Console.WriteLine($"Kodaira=[{tate.Select(e => e.kp).Glue(", ")}] Cp=[{tate.Select(e => e.cp).Glue(", ")}]");
    Console.WriteLine($"Conductor={N}");
    Console.WriteLine($"EllAn = [{ellAn.AscendingByKey().GlueMap(", ", "{0}:{1}")}]");
    Console.WriteLine($"L(E,1) = {LE1}");
    Console.WriteLine();

    return (N, ellAn, tate, Ell);
}

int EllConductor(BigInteger[] curve)
{
    var (a1, a2, a3, a4, a6) = curve.Select(i => new Rational(i)).Deconstruct();
    var E = new EllCoefs<Rational>(a1, a2, a3, a4, a6);
    E.Show();

    var dec = PrimesDec(E.Disc.Absolute.Num);
    var seq = dec.Keys.Select(p => TateAlgorithm(E, p)).ToArray();
    var N = seq.Select(e => e.p.Pow(e.fp)).Aggregate((pi, pj) => pi * pj);
    Console.WriteLine($"Kodaira=[{seq.Select(e => e.kp).Glue(", ")}] Cp=[{seq.Select(e => e.cp).Glue(", ")}]");
    Console.WriteLine($"Conductor={N}");
    Console.WriteLine();
    return N;
}

void testConductor()
{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    EllConductor([0, 0, 0, -1, 0]);
    EllConductor([0, 0, 0, 1, 0]);
    EllConductor([0, 0, 0, -5, 0]);
    EllConductor([0, 0, 0, 5, 0]);
    EllConductor([0, 0, 0, -25, 0]);
    EllConductor([0, 0, 0, -961, 0]);
    EllConductor([0, -1, 1, -10, -20]);
    EllConductor([0, 1, 0, 16, 180]);
    EllConductor([1, 1, 0, -22, -44]);
    EllConductor([1, -1, 1, -180, 1047]);
    EllConductor([0, 0, 0, -3, -18]);
    EllConductor([0, 0, 0, -123, -522]);
    EllConductor([0, -1, 1, 444, -826]);
    EllConductor([0, 1, 0, -5, 7]);
    EllConductor([1, -1, 0, -18, -81]);
    EllConductor([1, -1, 0, -17, 16]);
}

void testEllAn()
{
    var maxn = 32;
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;

    // Rank 0
    EllInfos([0, 0, 0, -1, 0], maxn);
    EllInfos([0, 0, 0, 1, 0], maxn);
    EllInfos([0, -1, 1, 0, 0], maxn);
    EllInfos([0, -1, 1, -10, -20], maxn);
    EllInfos([0, 0, 0, -3, -18], maxn);
    EllInfos([0, 0, 0, -123, -522], maxn);
    EllInfos([0, 1, 0, 16, 180], maxn);
}

double[] PrCoefs =
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

// RANK OF ELLIPTIC CURVES AND THE BIRCH
// SWINNERTON-DYER CONJECTURE
// J. HSU, S. MILLER
// Department of Mathematics
// Princeton University
double factorial(int k)
{
    var f = 1.0;
    for (int i = 1; i <= k; i++)
        f *= i;

    return f;
}

double Pr(int r, double x) => (r + 1).SeqLazy().Sum(k => PrCoefs[r - k] * double.Pow(x, k) / factorial(k));

double G(int r, double x)
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

void EllAnalyticRank(BigInteger[] curve)
{
    var (a1, a2, a3, a4, a6) = curve.Select(i => new Rational(i)).Deconstruct();
    var E = new EllCoefs<Rational>(a1, a2, a3, a4, a6);
    var Ell = new EllGroup<Rational>(a1, a2, a3, a4, a6);
    var dec = PrimesDec(E.Disc.Absolute.Num);
    var tate = dec.Keys.Select(p => TateAlgorithm(E, p)).ToArray();
    var N = tate.Select(e => new Rational(e.p).Pow(e.fp)).Aggregate((pi, pj) => pi * pj);
    var maxn = (int)double.Round(2 * double.Sqrt(N));
    var ellAn = EllAn(Ell, N, maxn);

    E.Show();
    Console.WriteLine($"Model {E.ModelStr}");
    Console.WriteLine($"Kodaira=[{tate.Select(e => e.kp).Glue(", ")}] Cp=[{tate.Select(e => e.cp).Glue(", ")}]");
    Console.WriteLine($"Conductor={N}");
    Console.WriteLine($"EllAn = [{ellAn.Where(e => e.Key <= 32).AscendingByKey().GlueMap(", ", "{0}:{1}")}]");

    var X = 2 * double.Pi / double.Sqrt(N);

    var r = 0;
    var L0 = maxn.SeqLazy(1).Sum(n => ellAn[n] * G(0, X * n) / n);
    var L1 = maxn.SeqLazy(1).Sum(n => ellAn[n] * G(1, X * n) / n);

    if (double.Abs(L0) > 0.001 && double.Abs(L1) > 0.001)
    {
        Console.WriteLine($"L^(0)(E, 1) = {L0}");
        Console.WriteLine($"L^(1)(E, 1) = {L1}");
        
        Console.WriteLine($"Analytic Rank = 0 or 1"); // TODO: sign of the functional equation
    }
    else
    {
        if (double.Abs(L0) < 0.001)
        {
            r = 2;
            Console.WriteLine($"L^(0)(E, 1) = {L0}");
        }
        else
        {
            r = 3;
            Console.WriteLine($"L^(1)(E, 1) = {L1}");
        }

        for (; r < 8; r += 2)
        {
            var Lr = factorial(r) * maxn.SeqLazy(1).Sum(n => ellAn[n] * G(r, X * n) / n);
            Console.WriteLine($"L^({r})(E, 1) = {Lr}");
            if (double.Abs(Lr) > 0.001)
                break;
        }
        
        Console.WriteLine($"Analytic Rank = {r}");
    }

    Console.WriteLine();
}

{
    Ring.DisplayPolynomial = MonomDisplay.StarCaret;
    EllAnalyticRank([0, 0, 0, -1, 0]);
    EllAnalyticRank([0, 0, 0, 1, 0]);
    EllAnalyticRank([0, -1, 1, 0, 0]);
    EllAnalyticRank([0, -1, 1, -10, -20]);
    EllAnalyticRank([0, 0, 0, -3, -18]);
    EllAnalyticRank([0, 0, 0, -123, -522]);
    EllAnalyticRank([0, 1, 0, 16, 180]);
    
    EllAnalyticRank([0, 0, 0, -5, 0]);
    EllAnalyticRank([0, 0, 0, -7, 0]);
    EllAnalyticRank([0, 0, 0, -25, 0]);
    EllAnalyticRank([0, 0, 0, -961, 0]);
    EllAnalyticRank([0, 0, 1, -1, 0]);
    
    EllAnalyticRank([0, 0, 0, -34 * 34, 0]);
    EllAnalyticRank([0, 0, 0, -41 * 41, 0]);
    EllAnalyticRank([0, 0, 0, -25, 25]);
    EllAnalyticRank([0, 0, 0, -81, 81]);
    EllAnalyticRank([1, 1, 0, -36, 36]);
    EllAnalyticRank([1, 1, 0, -49, 49]);
    EllAnalyticRank([0, 0, 1, -7, 6]);
    EllAnalyticRank([1, 1, 0, -87, 225]);
    EllAnalyticRank([1, 1, 0, -104, 276]);
    EllAnalyticRank([1, -1, 0, -79, 289]);
    EllAnalyticRank([0, 0, 1, -79, 342]);
}

public record TateAlgo(int p, int n, string kp, int fp, int cp);

public struct EllCoefs<K> where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
{
    public EllCoefs(K a1, K a2, K a3, K a4, K a6)
    {
        Model = (a1, a2, a3, a4, a6);

        var b2 = a1 * a1 + 4 * a2;
        var b4 = a1 * a3 + 2 * a4;
        var b6 = a3 * a3 + 4 * a6;
        var b8 = a1 * a1 * a6 - a1 * a3 * a4 + 4 * a2 * a6 + a2 * a3 * a3 - a4 * a4;
        B_Invariants = (b2, b4, b6, b8);

        var c4 = b2 * b2 - 24 * b4;
        var c6 = -b2.Pow(3) + 36 * b2 * b4 - 216 * b6;
        C_Invariants = (c4, c6);

        var disc = -b2 * b2 * b8 - 8 * b4.Pow(3) - 27 * b6 * b6 + 9 * b2 * b4 * b6;
        var j = c4.Pow(3) / disc;
        (J_Invariant, Disc) = (j, disc);
    }

    public (K a1, K a2, K a3, K a4, K a6) Model { get; }
    public (K b2, K b4, K b6, K b8) B_Invariants { get; }
    public (K c4, K c6) C_Invariants { get; }
    public K J_Invariant { get; }

    public (K a1, K a2, K a3, K a4, K a6, K b2, K b4, K b6, K b8, K c4, K c6, K j) ModelAndInvariants
    {
        get
        {
            var (a1, a2, a3, a4, a6) = Model;
            var (b2, b4, b6, b8) = B_Invariants;
            var (c4, c6) = C_Invariants;
            var j = J_Invariant;
            return (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, j);
        }
    }

    public K Disc { get; }

    public string Eq
    {
        get
        {
            var (x, y) = Ring.Polynomial(Disc, "x", "y").Deconstruct();
            var (a1, a2, a3, a4, a6) = Model;
            var lhs = y * y + a1 * x * y + a3 * y;
            var rhs = x.Pow(3) + a2 * x * x + a4 * x + a6;
            return $"Ellptic curve {lhs} = {rhs}";
        }
    }

    public string ModelStr => $"[{Model}]".Replace("(", "").Replace(")", "");

    public string B_InvariantsStr
    {
        get
        {
            var (b2, b4, b6, b8) = B_Invariants;
            return $"B Invariants b2={b2} b4={b4} b6={b6} b8={b8}";
        }
    }

    public string C_InvariantsStr
    {
        get
        {
            var (c4, c6) = C_Invariants;
            return $"C Invariants c4={c4} c6={c6}";
        }
    }

    public void Show()
    {
        Console.WriteLine(Eq);
        Console.WriteLine($"Disc = {Disc}");
        Console.WriteLine(B_InvariantsStr);
        Console.WriteLine(C_InvariantsStr);
        Console.WriteLine($"J Invariant j={J_Invariant}");
    }

    public EllCoefs<K> Transform(K r, K s, K t, K u)
    {
        var (a1, a2, a3, a4, a6) = Model;
        var a01 = (a1 + 2 * s) / u;
        var a02 = (a2 - s * a1 + 3 * r - s * s) / u.Pow(2);
        var a03 = (a3 + r * a1 + 2 * t) / u.Pow(3);
        var a04 = (a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1 + 3 * r * r - 2 * s * t) / u.Pow(4);
        var a06 = (a6 + r * a4 + r * r * a2 + r.Pow(3) - t * a3 - t * t - r * t * a1) / u.Pow(6);
        return new(a01, a02, a03, a04, a06);
    }

    public EllCoefs<K> Transform(int r, int s, int t, int u)
    {
        var o = Disc.One;
        return Transform(r * o, s * o, t * o, u * o);
    }
}