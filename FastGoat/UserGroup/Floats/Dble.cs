using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Floats;

public readonly struct Dble : IElt<Dble>, IRingElt<Dble>, IFieldElt<Dble>, IVsElt<Rational, Dble>, 
    IFloatElt<Dble>, IFixedPrecisionElt<Dble>
{
    public override bool Equals(object? obj)
    {
        return obj is Dble other && Equals(other);
    }

    public static double Eps => double.Pow(10, -Digits);
    public static Dble DbleZero() => new(0.0);
    public static Dble DbleOne() => new(1.0);
    public double K { get; }

    public Dble(double k)
    {
        K = double.Abs(k) < Eps * 10 ? 0 : k;
        Hash = K.GetHashCode();
    }

    public bool Equals(Dble other) => double.Abs(K - other.K) < Eps * 100;

    public int CompareTo(Dble other) => K.CompareTo(other.K);

    public int Hash { get; }
    public bool IsZero() => double.Abs(K) < Eps * 100;

    public Dble Zero => new(0);
    public Dble One => new(1);
    public Dble Add(Dble e) => new(K + e.K);

    public Dble Sub(Dble e) => new(K - e.K);

    public Dble Opp() => new(-K);

    public Dble Mul(Dble e) => new(K * e.K);

    public (Dble quo, Dble rem) Div(Dble e) => (new(K / e.K), Zero);

    public Dble Mul(int k) => new(K * k);

    public Dble Pow(int k)
    {
        if (k < 0)
            return Inv().Pow(-k);

        var d0 = 1.0;
        for (int i = 0; i < k; ++i)
            d0 *= K;

        return new(d0);
    }

    public static Dble operator +(Dble a, Dble b) => a.Add(b);

    public static Dble operator +(int a, Dble b) => new(a + b.K);

    public static Dble operator +(Dble a, int b) => new(a.K + b);

    public static Dble operator -(Dble a) => a.Opp();

    public static Dble operator -(Dble a, Dble b) => a.Sub(b);

    public static Dble operator -(int a, Dble b) => new(a - b.K);

    public static Dble operator -(Dble a, int b) => a.Sub(new(b));

    public static Dble operator *(Dble a, Dble b) => a.Mul(b);

    public static Dble operator *(int a, Dble b) => b.Mul(a);

    public static Dble operator *(Dble a, int b) => a.Mul(b);

    public static Dble operator /(Dble a, Dble b) => a.Div(b).quo;

    public static Dble operator /(Dble a, int b) => new(a.K / b);

    public int P => 0;
    public Dble Inv() => new(1.0 / K);
    public int Sign => double.Sign(K);

    public bool Invertible() => !IsZero();
    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public Dble Absolute => new(Sign * K);
    public Dble Sqrt() => new(double.Sqrt(K));
    public static Dble From<T>(T e) where T : IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
    {
        if (e is BigReal e0)
            return e0.ToDble;

        if (e is Dcml e1)
            return new((double)e1.K);

        if (e is Rational e2)
            return new(e2.ToDouble);

        if (e is Dble e3)
            return e3;

        throw new ArgumentException();
    }

    public static int Digits => 17;

    public Dble RoundEven => Round(this, 0);
    public BigReal ToBigReal(int O) => BigReal.FromString($"{K:E18}", O);

    public Dble KMul(Rational k) => new(K * k.ToDouble);

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{K}";

    public static Dble NthRoot(Dble a, int n) => new(double.RootN(a.K, n));
    public static Dble Sqrt(Dble a) => new(double.Sqrt(a.K));

    public static Dble operator /(int a, Dble b) => a * b.Inv();

    public static double Abs(Dble t) => (double)double.Abs(t.K);
    public static Dble Round(Dble e, int d, MidpointRounding mode = MidpointRounding.ToEven) 
        => new(double.Round(e.K, d, mode));

    public static bool IsValuedField => true;

    public static implicit operator double(Dble e) => e.K;

    public static implicit operator Dble(Rational e) => new((double)e.ToDouble);

    public static Dble operator +(Dble a, Rational b) => a.Add(b);

    public static Dble operator +(Rational a, Dble b) => b.Add(a);

    public static Dble operator -(Dble a, Rational b) => a.Sub(b);

    public static Dble operator -(Rational a, Dble b) => new Dble((double)a.ToDouble).Sub(b);

    public static Dble operator *(Dble a, Rational b) => a.Mul(new Dble((double)b.ToDouble));

    public static Dble operator *(Rational a, Dble b) => new Dble((double)a.ToDouble).Mul(b);

    public static Dble operator /(Dble a, Rational b) => a.Div(b).quo;

    public static bool operator ==(Dble a, Dble b) => a.Equals(b);

    public static bool operator !=(Dble a, Dble b) => !(a == b);
    public static bool operator <(Dble a, Dble b) => a.K < b.K;

    public static bool operator >(Dble a, Dble b) => a.K > b.K;
    public static bool operator <=(Dble a, Dble b) => (a < b) || (a == b);

    public static bool operator >=(Dble a, Dble b) => (a > b) || (a == b);
    public static Dble Min(Dble a, Dble b) => a <= b ? a : b;
    public static Dble Max(Dble a, Dble b) => a < b ? b : a;
}