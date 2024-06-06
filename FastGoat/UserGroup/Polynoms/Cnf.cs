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
        CnfOne = new Cnf(1).One;
    }

    public static char RootsOfUnit { get; set; } = 'ξ';

    public static Cnf Nth(int k) => new(k);

    public static Cnf CnfZero { get; }
    public static Cnf CnfOne { get; }
    public static Cnf I { get; }
    public static KPoly<Cnf> X => FG.KPoly('X', CnfOne);
    public int N { get; }
    public EPoly<Rational> E { get; }

    public bool IsInteger => E.Degree == 0 && E[0].IsInteger();
    public bool IsPositiveInteger => IsInteger && E[0].Sign == 1;

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

    public double Phase
    {
        get
        {
            var c = Complex.FromPolarCoordinates(1, 2 * Double.Pi / N);
            var c0 = E.Poly.Substitute(c);
            return c0.Phase;
        }
    }

    public Cnf(int n)
    {
        E = FG.CyclotomicEPoly(n);
        N = n;
        Hash = 1;
    }

    public Cnf(int n, EPoly<Rational> e)
    {
        E = e;
        N = n;
        Hash = 1;
    }

    public Complex ToComplex => Complex.FromPolarCoordinates(Module, Phase);

    public int GetHashCodeSlow() => E.GetHashCodeSlow();
    public bool Equals(Cnf other) => (this - other).IsZero();

    public int CompareTo(Cnf other)
    {
        var ce = ToComplex;
        var co = other.ToComplex;
        var compMod = ce.Magnitude.CompareTo(co.Magnitude);
        if (compMod != 0)
            return compMod;

        return ce.Phase.CompareTo(co.Phase);
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
    public bool Invertible() => true;

    public (Cnf quo, Cnf rem) Div(Cnf e) => (Mul(e.Inv()), Zero);

    public Cnf Mul(int k) => new(N, k * E);
    public Cnf Mul(Rational k) => new(N, k * E);

    public Cnf Pow(int k) => new(N, E.Pow(k));
    public Cnf Simplify() => FG.CnfBasis(N).Simplify(this).Item1;

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var (cf, p0) = FG.CnfBasis(N).Simplify(this);
        if (N == 4)
            p0 = p0.Div(cf.E.F).rem;
        
        var letter = cf.N == 4 ? "I" : $"{RootsOfUnit}{cf.N}";
        var ind = Ring.Indeterminates(letter);
        ind.SetOrder(MonomOrder.RevLex);

        var polStr = p0.ToPolynomial(ind, ind[0]);
        return $"{polStr}";
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
        return c0.Magnitude;
    }

    public static bool IsValuedField => true;
}