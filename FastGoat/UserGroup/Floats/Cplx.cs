using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Floats;

public readonly struct Cplx : IElt<Cplx>, IRingElt<Cplx>, IFieldElt<Cplx>, IVsElt<Rational, Cplx>, 
    IFloatElt<Cplx>, IFixedPrecisionElt<Cplx>
{

    public static bool IsComplex => true;
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
        var a0 = double.Abs(k.Real) < Eps ? 0.0 : k.Real;
        var b0 = double.Abs(k.Imaginary) < Eps ? 0.0 : k.Imaginary;
        var k0 = a0 + Complex.ImaginaryOne * b0;
        var nInf = Double.Abs(Double.Max(Double.Abs(k0.Real), Double.Abs(k0.Imaginary)));
        K = nInf < Eps ? Complex.Zero : k0;
        Hash = K.GetHashCode();
    }

    public bool Equals(Cplx other) => Sub(other).NormInf < Eps * 10000;
    public override bool Equals(object? obj)
    {
        return obj is Cplx other && Equals(other);
    }

    public int CompareTo(Cplx other)
    {
        var mag = K.Magnitude.CompareTo(other.K.Magnitude);
        if (Double.Abs(K.Magnitude - other.K.Magnitude) > Eps)
            return mag;

        return Phase.CompareTo(other.Phase);
    }

    public int Hash { get; init; }
    public Complex ToComplex => RealPart + Complex.ImaginaryOne * ImaginaryPart;
    public bool IsZero() => K.Magnitude < Eps * 10000;
    public bool IsReal() => Double.Abs(K.Imaginary) < Eps * 10000;
    public bool IsImaginary() => Double.Abs(K.Real) < Eps * 10000;
    public BigCplx ToBigComplex(int O = 40) => new(BigReal.FromDouble(RealPart, O), BigReal.FromDouble(ImaginaryPart, O));

    public static Cplx I => new(Complex.ImaginaryOne);

    public static Cplx CZero => new(Complex.Zero);

    public static Cplx COne => new(Complex.One);
    public Cplx RoundEven => new(double.Round(RealPart) + Complex.ImaginaryOne * double.Round(ImaginaryPart));
    public Cplx Absolute => new(Magnitude);
    public Cplx Absolute2 => new(RealPart * RealPart + ImaginaryPart * ImaginaryPart);
    public Cplx Conj => new(Complex.Conjugate(K));

    public int Sign => K.Real == 0.0 ? 1 : double.Sign(K.Real);
    public double ToDouble =>  throw new NotImplementedException();
    public static Cplx Pi(int o = 50) => new(Double.Pi);

    public static Cplx Sqrt(Cplx e) => new(Complex.Sqrt(e.K));
    public static Cplx NthRoot(Cplx r, int n) => new(Complex.Pow(r.ToComplex, 1.0 / n));

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
        var epsSingle = double.Pow(10, -Digits / 2.0);
        if (K.Magnitude < epsSingle)
            return "0";
        if ((K - Complex.One).Magnitude < epsSingle)
            return "1";
        if ((K + Complex.One).Magnitude < epsSingle)
            return "-1";

        var sep = (Ring.DisplayPolynomial & MonomDisplay.Star) == MonomDisplay.Star ? "*" : "·";

        var a0 = $"{Double.Round(K.Real, Digits / 2)}";
        var b0 = Double.Abs(K.Imaginary - 1.0) < epsSingle
            ? "I"
            : Double.Abs(K.Imaginary + 1.0) < epsSingle
                ? "-I"
                : $"{Double.Round(K.Imaginary, Digits / 2)}{sep}I";

        if (Double.Abs(K.Imaginary) < epsSingle)
            return a0;

        if (Double.Abs(K.Real) < epsSingle)
            return b0;

        return $"({a0} + {b0})".Replace("+ -", "- ");
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
    public static Cplx From<T>(T e) where T : IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
    {
        if (e is BigReal e0)
            return new(e0.ToDouble);
        
        if (e is Dcml e1)
            return new((double)e1.K);

        if (e is Rational e2)
            return new(e2.ToDouble);

        if (e is Dble e3)
            return new(e3.K);

        if (e is BigCplx e4)
            return new(e4.ToComplex);

        throw new ArgumentException();
    }

    public static int Digits => 17;
    public static double Eps => double.Pow(10, -Digits);
    public static bool operator ==(Cplx a, Cplx b) => a.Equals(b);

    public static bool operator !=(Cplx a, Cplx b) => !a.Equals(b);

    public static bool operator <(Cplx a, Cplx b) => a.CompareTo(b) == -1;

    public static bool operator >(Cplx a, Cplx b) =>  a.CompareTo(b) == 1;

    public static bool operator <=(Cplx a, Cplx b) =>  a.CompareTo(b) <= 0;

    public static bool operator >=(Cplx a, Cplx b) => a.CompareTo(b) >= 0;

    public static Cplx Min(Cplx a, Cplx b) => a.CompareTo(b) <= 0 ? a : b;

    public static Cplx Max(Cplx a, Cplx b)=> a.CompareTo(b) >= 0 ? a : b;
}