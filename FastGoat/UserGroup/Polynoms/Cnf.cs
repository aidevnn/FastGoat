using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public struct Cnf : IElt<Cnf>, IRingElt<Cnf>, IFieldElt<Cnf>
{
    static Cnf()
    {
        I = new Cnf(4);
        CnfZero = new Cnf(1).Zero;
        CnfOne = new Cnf(1);
    }

    public static Cnf Nth(int k) => new(k);

    public static Cnf CnfZero { get; }
    public static Cnf CnfOne { get; }
    public static Cnf I { get; }
    public int N { get; }
    public EPoly<Rational> E { get; }

    public Cnf Re
    {
        get
        {
            var a = E.X;
            var e = (E + E.Substitute(a.Inv())) / 2;
            return new(N, e);
        }
    }

    public Cnf Im => (this - Re) / I;
    public Cnf Conj => new(N, E.Substitute(E.X.Inv()));

    public Cnf Module2
    {
        get
        {
            var re = Re;
            var im = Im;
            return re * re + im * im;
        }
    }

    public double Module => Abs(this);

    public Cnf(int n)
    {
        E = FG.CyclotomicEPoly(n);
        N = n;
        Hash = (N, E.Hash).GetHashCode();
    }

    private Cnf(int n, EPoly<Rational> e)
    {
        E = e;
        N = n;
        Hash = (N, E.Hash).GetHashCode();
    }

    public bool Equals(Cnf other) => (this - other).IsZero();

    public int CompareTo(Cnf other)
    {
        var e = (this - other).E;
        return e.CompareTo(e.Zero);
    }

    public int Hash { get; }
    public bool IsZero() => E.IsZero();

    public Cnf Zero => new(1, FG.CyclotomicEPoly(1).Zero);
    public Cnf One => new(1, FG.CyclotomicEPoly(1).One);

    public int P => 0;

    public Cnf Add(Cnf e)
    {
        var lcm = (N * e.N) / (IntExt.Gcd(N, e.N));
        var a = FG.CyclotomicEPoly(lcm);
        var a0 = a.Pow(lcm / N);
        var a1 = a.Pow(lcm / e.N);
        return new(lcm, E.Substitute(a0) + e.E.Substitute(a1));
    }

    public Cnf Opp() => new(N, -E);

    public Cnf Sub(Cnf e) => Add(e.Opp());

    public Cnf Mul(Cnf e)
    {
        var lcm = (N * e.N) / (IntExt.Gcd(N, e.N));
        var a = FG.CyclotomicEPoly(lcm);
        var a0 = a.Pow(lcm / N);
        var a1 = a.Pow(lcm / e.N);
        return new(lcm, E.Substitute(a0) * e.E.Substitute(a1));
    }

    public Cnf Inv() => new(N, E.Inv());

    public (Cnf quo, Cnf rem) Div(Cnf e) => (Mul(e.Inv()), Zero);

    public Cnf Mul(int k) => new(N, k * E);
    public Cnf Mul(Rational k) => new(N, k * E);

    public Cnf Pow(int k) => new(N, E.Pow(k));
    public Cnf Simplify() => Simplify(this);

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var cf = Simplify(this);
        var letter = cf.N == 4 ? "I" : $"ξ{cf.N}";
        var ind = Ring.Indeterminates(letter);
        var poly = cf.E.Poly.ToPolynomial(ind, ind[0]);
        return $"{poly}";
    }

    public static Cnf Simplify(Cnf c)
    {
        if (c.IsZero())
            return CnfZero;

        var a = c.E.X;
        var n0 = c.N;
        var a0 = a.One;
        for (int i = 0; i < n0; i++)
        {
            var cfe = c.E;
            var t = cfe.Div(a0).quo;
            if (t.Poly.Degree == 0)
            {
                var gcd = IntExt.Gcd(n0, i);
                var n1 = n0 / gcd;
                var i1 = i / gcd;
                var a1 = FG.CyclotomicEPoly(n1).Pow(i1);
                var a2 = new Cnf(n1, t.Poly[0] * a1);
                return a2;
            }

            a0 *= a;
        }

        return c;
    }

    public static Cnf operator +(Cnf a, Cnf b) => a.Add(b);

    public static Cnf operator +(int a, Cnf b) => b.Add(b.One.Mul(a));

    public static Cnf operator +(Cnf a, int b) => a.Add(a.One.Mul(b));

    public static Cnf operator +(Rational a, Cnf b) => b.Add(b.One.Mul(a));

    public static Cnf operator +(Cnf a, Rational b) => a.Add(a.One.Mul(b));

    public static Cnf operator -(Cnf a) => a.Opp();

    public static Cnf operator -(Cnf a, Cnf b) => a.Sub(b);

    public static Cnf operator -(int a, Cnf b) => b.One.Mul(a).Sub(b);

    public static Cnf operator -(Cnf a, int b) => a.Sub(a.One * b);

    public static Cnf operator -(Rational a, Cnf b) => b.One.Mul(a).Sub(b);

    public static Cnf operator -(Cnf a, Rational b) => a.Sub(a.One * b);

    public static Cnf operator *(Cnf a, Cnf b) => a.Mul(b);

    public static Cnf operator *(int a, Cnf b) => b.Mul(a);

    public static Cnf operator *(Cnf a, int b) => a.Mul(b);

    public static Cnf operator *(Rational a, Cnf b) => b.Mul(a);

    public static Cnf operator *(Cnf a, Rational b) => a.Mul(b);

    public static Cnf operator /(Cnf a, Cnf b) => a.Div(b).quo;

    public static Cnf operator /(Cnf a, int b) => a.Div(a.One * b).quo;

    public static Cnf operator /(int a, Cnf b) => b.Inv().Mul(a);

    public static Cnf operator /(Cnf a, Rational b) => a.Div(a.One * b).quo;

    public static Cnf operator /(Rational a, Cnf b) => b.Inv().Mul(a);

    public static double Abs(Cnf t)
    {
        var c = Complex.FromPolarCoordinates(1, 2 * Double.Pi / t.N);
        var c0 = t.E.Poly.Substitute(c);
        return Complex.Abs(c0);
    }

    public static bool IsValuedField => true;
}