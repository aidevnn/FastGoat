using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Integers;
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

(string Kp, int fp, int cp) TateAlgorithm(EllCurveCoefs<Rational> E, int p)
{
    var n = Ord(p, E.Disc.Num);
    if (n == 0)
        return ("I0", 0, 1);

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

        return ($"I{n}", 1, cp);
    }

    if (!a6.Mod(p * p).IsZero())
        return ("II", n, 1);

    if (!b8.Mod(p * p * p).IsZero())
        return ("III", n - 1, 2);

    if (!b6.Mod(p * p * p).IsZero())
    {
        if (QuadRoot(a1.One, a3 / p, -a6 / (p * p), p))
            cp = 3;
        else
            cp = 1;

        return ("IV", n - 2, cp);
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
        return ("I*0", n - 4, 1 + NCubicRoots(b, c, d, p));
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
            // Console.WriteLine(new { dbg = 1, p, cp, m, mx, my, xa2, xa3, xa4, xa6 });
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
                // Console.WriteLine(new { dbg = 2, p, cp, m, mx, my, xa2, xa3, xa4, xa6 });
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

        return ($"I*{m}", n - m - 4, cp);
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

            return ("IV*", n - 6, cp);
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
                return ("III*", n - 7, 2);
            else if (!a6.Mod(p.Pow(6)).IsZero())
                return ("II*", n - 8, 1);
            else
            {
                Console.WriteLine("restart");
                return TateAlgorithm(Etmp.Transform(0, 0, 0, p), p);
            }
        }
    }
}

int EllConductor(BigInteger[] curve)
{
    var (a1, a2, a3, a4, a6) = curve.Select(i => new Rational(i)).Deconstruct();
    var E = new EllCurveCoefs<Rational>(a1, a2, a3, a4, a6);
    E.Show();

    var dec = PrimesDec(E.Disc.Absolute.Num);
    var N = 1;
    var seqKp = new List<string>();
    var seqCp = new List<int>();
    foreach (var p in dec.Keys)
    {
        var (kp, fp, cp) = TateAlgorithm(E, p);
        seqKp.Add(kp);
        seqCp.Add(cp);
        N *= p.Pow(fp);
    }

    Console.WriteLine($"Kodaira=[{seqKp.Glue(", ")}] Cp=[{seqCp.Glue(", ")}]");
    Console.WriteLine($"Conductor={N}");
    Console.WriteLine();
    return N;
}

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

public struct EllCurveCoefs<K> where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
{
    public EllCurveCoefs(K a1, K a2, K a3, K a4, K a6)
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
        var (x, y) = Ring.Polynomial(Disc, "x", "y").Deconstruct();
        var (a1, a2, a3, a4, a6) = Model;
        var lhs = y * y + a1 * x * y + a3 * y;
        var rhs = x.Pow(3) + a2 * x * x + a4 * x + a6;
        Console.WriteLine($"Ellptic curve {lhs} = {rhs}");
        Console.WriteLine($"Disc = {Disc}");
        Console.WriteLine(B_InvariantsStr);
        Console.WriteLine(C_InvariantsStr);
        Console.WriteLine($"J Invariant j={J_Invariant}");
    }

    public EllCurveCoefs<K> Transform(K r, K s, K t, K u)
    {
        var (a1, a2, a3, a4, a6) = Model;
        var a01 = (a1 + 2 * s) / u;
        var a02 = (a2 - s * a1 + 3 * r - s * s) / u.Pow(2);
        var a03 = (a3 + r * a1 + 2 * t) / u.Pow(3);
        var a04 = (a4 - s * a3 + 2 * r * a2 - (t + r * s) * a1 + 3 * r * r - 2 * s * t) / u.Pow(4);
        var a06 = (a6 + r * a4 + r * r * a2 + r.Pow(3) - t * a3 - t * t - r * t * a1) / u.Pow(6);
        return new(a01, a02, a03, a04, a06);
    }

    public EllCurveCoefs<K> Transform(int r, int s, int t, int u)
    {
        var o = Disc.One;
        return Transform(r * o, s * o, t * o, u * o);
    }
}