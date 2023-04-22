using System.Net;
using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Polynoms;

public readonly struct Cplx : IElt<Cplx>, IRingElt<Cplx>, IFieldElt<Cplx>, IVsElt<Rational, Cplx>
{
    public const int Digits = 8;
    public const double EpsSingle = 1e-8;
    public const double EpsDouble = 1e-15;
    public Complex K { get; }
    public double RealPart => K.Real;
    public double ImaginaryPart => K.Imaginary;
    public double Magnitude => K.Magnitude;

    public double Phase
    {
        get
        {
            var p = K.Phase;
            return p > 0 ? p : Double.Pi * 2 + p;
        }
    }
    public double NormInf => Double.Max(Double.Abs(K.Real), Double.Abs(K.Imaginary));
    public double Norm1 => Double.Abs(K.Real) + Double.Abs(K.Imaginary);
    public Cplx(Complex k)
    {
        var a0 = double.Abs(k.Real) < EpsDouble ? 0.0 : k.Real;
        var b0 = double.Abs(k.Imaginary) < EpsDouble ? 0.0 : k.Imaginary;
        var k0 = a0 + Complex.ImaginaryOne * b0;
        var nInf = Double.Abs(Double.Max(Double.Abs(k0.Real), Double.Abs(k0.Imaginary)));
        K = nInf < EpsDouble ? Complex.Zero : k0;
        Hash = K.GetHashCode();
    }

    public bool Equals(Cplx other) => Sub(other).NormInf < EpsSingle;

    public int CompareTo(Cplx other)
    {
        var mag = K.Magnitude.CompareTo(other.K.Magnitude);
        if (Double.Abs(K.Magnitude - other.K.Magnitude) > EpsSingle)
            return mag;

        return Phase.CompareTo(other.Phase);
    }

    public int Hash { get; init; }
    public bool IsZero() => K.Magnitude < EpsSingle;
    public bool IsReal() => Double.Abs(K.Imaginary) < EpsSingle;
    public bool IsImaginary() => Double.Abs(K.Real) < EpsSingle;

    public static Cplx I => new(Complex.ImaginaryOne);

    public static Cplx CZero => new(Complex.Zero);

    public static Cplx COne => new(Complex.One);
    public static Cplx Sqrt(Cplx e) => new(Complex.Sqrt(e.K));
    public static KPoly<Cplx> X => new('X', CZero, new[] { CZero, COne });
    public Cplx Zero => new(Complex.Zero);
    public Cplx One => new(Complex.One);
    public Cplx Add(Cplx e) => new(K + e.K);

    public Cplx Sub(Cplx e) => new(K - e.K);

    public Cplx Opp() => new(-K);

    public Cplx Mul(Cplx e) => new(K * e.K);

    public (Cplx quo, Cplx rem) Div(Cplx e) => (new(K / e.K), Zero);

    public Cplx Mul(int k) => new(K * k);

    public Cplx Pow(int k) => new(Complex.Pow(K, k));

    public static implicit operator Cplx(Rational r) => new(r + 0.0);

    public static Cplx operator +(Cplx a, Cplx b) => a.Add(b);

    public static Cplx operator +(int a, Cplx b) => new(a + b.K);

    public static Cplx operator +(Cplx a, int b) => new(a.K + b);

    public static Cplx operator -(Cplx a) => a.Opp();

    public static Cplx operator -(Cplx a, Cplx b) => a.Sub(b);

    public static Cplx operator -(int a, Cplx b) => new(a - b.K);

    public static Cplx operator -(Cplx a, int b) => new(a.K - b);

    public static Cplx operator *(Cplx a, Cplx b) => a.Mul(b);

    public static Cplx operator *(int a, Cplx b) => b.Mul(a);

    public static Cplx operator *(Cplx a, int b) => a.Mul(b);

    public static Cplx operator /(Cplx a, Cplx b) => a.Div(b).quo;

    public static Cplx operator /(Cplx a, int b) => new(a.K / b);

    public int P => 0;
    public Cplx Inv() => new(1.0 / K);

    public bool Invertible() => IsZero();

    public static Cplx operator /(int a, Cplx b) => new(a / b.K);

    public static double Abs(Cplx t) => t.K.Magnitude;

    public static bool IsValuedField => true;

    public override int GetHashCode() => IsZero() ? 0 : IsReal() ? 1 : 2;

    public override string ToString()
    {
        if (K.Magnitude < EpsSingle)
            return "0";
        if ((K - Complex.One).Magnitude < EpsSingle)
            return "1";
        if ((K + Complex.One).Magnitude < EpsSingle)
            return "-1";

        var sep = (Ring.DisplayPolynomial & MonomDisplay.Star) == MonomDisplay.Star ? "*" : "·";

        var a0 = $"{Double.Round(K.Real, Digits)}";
        var b0 = Double.Abs(K.Imaginary - 1.0) < EpsSingle
            ? "I"
            : Double.Abs(K.Imaginary + 1.0) < EpsSingle
                ? "-I"
                : $"{Double.Round(K.Imaginary, Digits)}{sep}I";

        if (Double.Abs(K.Imaginary) < EpsSingle)
            return a0;

        if (Double.Abs(K.Real) < EpsSingle)
            return b0;

        return $"({a0} + {b0})";
    }

    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public Cplx KMul(Rational k) => ((Cplx)k) * this;

    public static Cplx operator +(Cplx a, Rational b) => a + ((Cplx)b);

    public static Cplx operator +(Rational a, Cplx b) => ((Cplx)a) + b;

    public static Cplx operator -(Cplx a, Rational b) => a - ((Cplx)b);

    public static Cplx operator -(Rational a, Cplx b) => ((Cplx)a) - b;

    public static Cplx operator *(Cplx a, Rational b) => a.KMul(b);

    public static Cplx operator *(Rational a, Cplx b) => b.KMul(a);

    public static Cplx operator /(Cplx a, Rational b) => a / ((Cplx)b);
}