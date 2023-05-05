using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public readonly struct BigCplx : IElt<BigCplx>, IRingElt<BigCplx>, IFieldElt<BigCplx>, IVsElt<Rational, BigCplx>
{
    public static bool IsValuedField => true;

    public static double Abs(BigCplx t) => t.ToComplex.Magnitude;
    public static BigCplx BcZero(int o = 40) => new(o);
    public static BigCplx BcOne(int o = 40) => BcZero(o).One;
    public static BigCplx BgI(int o = 40) => BcZero(o).I;

    public BigReal RealPart { get; }
    public BigReal ImaginaryPart { get; }
    public int O => RealPart.O;

    public BigCplx(int o)
    {
        RealPart = BigReal.BrZero(o);
        ImaginaryPart = BigReal.BrZero(o);
        Hash = (RealPart, ImaginaryPart).GetHashCode();
    }

    public BigCplx(BigReal re, BigReal im)
    {
        if (re.O != im.O)
            throw new($"Real part and Imaginary part must have the same digits precision");

        RealPart = re;
        ImaginaryPart = im;
        Hash = (RealPart, ImaginaryPart).GetHashCode();
    }

    public BigCplx Conj => new(RealPart, -ImaginaryPart);
    public BigReal NormInf => BigReal.Max(RealPart.Absolute, ImaginaryPart.Absolute);

    public bool Equals(BigCplx other)
    {
        var diffRe = (RealPart - other.RealPart);
        var diffIm = (ImaginaryPart - other.ImaginaryPart);
        return diffRe.IsZero() && diffIm.IsZero();
    }

    public int CompareTo(BigCplx other)
    {
        var mag = Magnitude2.CompareTo(other.Magnitude2);
        if (mag != 0)
            return mag;

        return Phase.CompareTo(other.Phase);
    }

    public int Hash { get; }

    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();

    public Complex ToComplex => new(RealPart.ToDouble, ImaginaryPart.ToDouble);

    public BigReal Magnitude2 => RealPart.Pow(2) + ImaginaryPart.Pow(2);

    public double Magnitude => ToComplex.Magnitude;
    public double Phase => ToComplex.Phase;

    public int P => 0;
    public bool Invertible() => !IsZero();
    public bool IsZero() => RealPart.IsZero() && ImaginaryPart.IsZero();
    public bool IsZero3d() => RealPart.IsZero3d() && ImaginaryPart.IsZero3d();
    public bool IsZero4d() => RealPart.IsZero4d() && ImaginaryPart.IsZero4d();

    public BigCplx Zero => new(BigReal.BrZero(O), BigReal.BrZero(O));
    public BigCplx One => new(BigReal.BrOne(O), BigReal.BrZero(O));
    public BigCplx I => new(BigReal.BrZero(O), BigReal.BrOne(O));
    public BigCplx Add(BigCplx e) => new(RealPart + e.RealPart, ImaginaryPart + e.ImaginaryPart);

    public BigCplx Sub(BigCplx e) => new(RealPart - e.RealPart, ImaginaryPart - e.ImaginaryPart);

    public BigCplx Opp() => new(RealPart.Opp(), ImaginaryPart.Opp());

    public BigCplx Mul(BigCplx e) => new(RealPart * e.RealPart - ImaginaryPart * e.ImaginaryPart,
        RealPart * e.ImaginaryPart + e.RealPart * ImaginaryPart);

    public (BigCplx quo, BigCplx rem) Div(BigCplx e) => (Mul(e.Inv()), Zero);

    public BigCplx Inv()
    {
        if (IsZero())
            throw new DivideByZeroException();

        var abs = Magnitude2;
        return new(RealPart / abs, -ImaginaryPart / abs);
    }

    public BigCplx Mul(int k) => new(RealPart * k, ImaginaryPart * k);

    public BigCplx KMul(Rational k)
    {
        var k0 = BigReal.FromRational(k, O);
        return new(RealPart * k0, ImaginaryPart * k0);
    }

    public BigCplx Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var e = this;
        return Enumerable.Repeat(e, k).Aggregate(One, (acc, f) => acc * f);
    }

    public BigCplx Pow10(int n) => One.Mul(10).Pow(n);
    public BigCplx ToBigCplx(int o) => new(RealPart.ToBigReal(o), ImaginaryPart.ToBigReal(o));

    public string ToSciForm()
    {
        var fmt = $"[{{0,-{O + 7}}} ; {{1,-{O + 7}}}]";
        return String.Format(fmt, RealPart.ToSciForm(), ImaginaryPart.ToSciForm());
    }

    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        if (IsZero())
            return "0";
        if ((this - 1).IsZero())
            return "1";
        if ((this + 1).IsZero())
            return "-1";

        var sep = (Ring.DisplayPolynomial & MonomDisplay.Star) == MonomDisplay.Star ? "*" : "·";

        var a0 = $"{RealPart.ToDouble:G15}";
        var b0 = (ImaginaryPart-1).IsZero()
            ? "I"
            : (ImaginaryPart+1).IsZero()
                ? "-I"
                : $"{ImaginaryPart.ToDouble:G15}{sep}I";

        if (ImaginaryPart.IsZero())
            return a0;

        if (RealPart.IsZero())
            return b0;

        return $"({a0} + {b0})";
    }

    public static BigCplx Round(BigCplx c, int digits) => new(BigReal.Round(c.RealPart, digits), BigReal.Round(c.ImaginaryPart, digits));

    public static BigCplx FromRational(Rational r, int O)
    {
        var re = BigReal.FromRational(r, O);
        return new(re, re.Zero);
    }

    public static BigCplx FromRational(Rational re, Rational im, int O) =>
        FromBigReal(BigReal.FromRational(re, O), BigReal.FromRational(im, O));
    public static BigCplx FromBigReal(BigReal re) => new(re, re.Zero);
    public static BigCplx FromBigReal(BigReal re, BigReal im) => new(re, im);

    public static BigReal MagnitudeBigReal(BigCplx r) => BigReal.Sqrt(r.Magnitude2);

    public static BigCplx Sqrt(BigCplx r) => NthRoot(r, 2);

    public static BigCplx NthRoot(BigCplx r, int n)
    {
        if (n == 0)
            return r.One;

        if (n < 0)
            return NthRoot(r, -n).Inv();

        var ai = r.Zero;
        var aj = new BigCplx(BigReal.FromDouble(Double.Pi, r.O), BigReal.FromDouble(Double.E, r.O));
        while (!(ai - aj).IsZero())
        {
            ai = aj;
            var aiPow = ai.Pow(n - 1);
            var num = aiPow * ai - r;
            var denom = n * aiPow;
            aj = ai - num / denom; // Newton iteration
        }

        return aj;
    }

    public static BigCplx operator +(BigCplx a, BigCplx b) => a.Add(b);

    public static BigCplx operator +(int a, BigCplx b) => b.One.Mul(a).Add(b);

    public static BigCplx operator +(BigCplx a, int b) => a.Add(a.One.Mul(b));

    public static BigCplx operator -(BigCplx a) => a.Opp();

    public static BigCplx operator -(BigCplx a, BigCplx b) => a.Sub(b);

    public static BigCplx operator -(int a, BigCplx b) => b.One.Mul(a).Sub(b);

    public static BigCplx operator -(BigCplx a, int b) => a.Sub(a.One.Mul(b));

    public static BigCplx operator *(BigCplx a, BigCplx b) => a.Mul(b);

    public static BigCplx operator *(int a, BigCplx b) => b.Mul(a);

    public static BigCplx operator *(BigCplx a, int b) => a.Mul(b);

    public static BigCplx operator /(BigCplx a, BigCplx b) => a.Div(b).quo;

    public static BigCplx operator /(BigCplx a, int b) => a.Div(a.One.Mul(b)).quo;

    public static BigCplx operator /(int a, BigCplx b) => b.Inv().Mul(a);

    public static BigCplx operator +(BigCplx a, Rational b) => new(a.RealPart + BigReal.FromRational(b, a.O), a.ImaginaryPart);

    public static BigCplx operator +(Rational a, BigCplx b) => new(b.RealPart + BigReal.FromRational(a, b.O), b.ImaginaryPart);

    public static BigCplx operator -(BigCplx a, Rational b) => new(a.RealPart - BigReal.FromRational(b, a.O), a.ImaginaryPart);

    public static BigCplx operator -(Rational a, BigCplx b) => new(b.RealPart - BigReal.FromRational(a, b.O), b.ImaginaryPart);

    public static BigCplx operator *(BigCplx a, Rational b)
    {
        var b0 = BigReal.FromRational(b, a.O);
        return new(a.RealPart * b0, a.ImaginaryPart * b0);
    }

    public static BigCplx operator *(Rational a, BigCplx b) 
    {
        var a0 = BigReal.FromRational(a, b.O);
        return new(b.RealPart * a0, b.ImaginaryPart * a0);
    }

    public static BigCplx operator /(BigCplx a, Rational b) 
    {
        var b0 = BigReal.FromRational(b, a.O);
        return new(a.RealPart / b0, a.ImaginaryPart / b0);
    }
}