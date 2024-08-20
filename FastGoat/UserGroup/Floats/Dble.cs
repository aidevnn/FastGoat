using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Floats;

public readonly struct Dble : IElt<Dble>, IRingElt<Dble>, IFieldElt<Dble>, IVsElt<Rational, Dble>
{
    public static double EpsDouble = 1e-14;
    public static Dble DbleZero() => new(0.0);
    public static Dble DbleOne() => new(1.0);
    public Double K { get; }

    public Dble(Double k)
    {
        K = Double.Abs(k) < EpsDouble ? 0 : k;
        Hash = K.GetHashCode();
    }

    public bool Equals(Dble other) => double.Abs(K - other.K) < EpsDouble;

    public int CompareTo(Dble other) => K.CompareTo(other.K);

    public int Hash { get; }
    public bool IsZero() => double.Abs(K) < EpsDouble;

    public Dble Zero => new(0);
    public Dble One => new(1);
    public Dble Add(Dble e) => new(K + e.K);

    public Dble Sub(Dble e) => new(K - e.K);

    public Dble Opp() => new(-K);

    public Dble Mul(Dble e) => new(K * e.K);

    public (Dble quo, Dble rem) Div(Dble e) => (new(K / e.K), Zero);

    public Dble Mul(int k) => new(K * k);

    public Dble Pow(int k) => new(Double.Pow(K, k));

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
    public Dble Inv() => new(1 / K);
    public int Sign => double.Sign(K);

    public bool Invertible() => !IsZero();
    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public Dble Absolute => new(Sign * K);
    public Dble RoundEven => Round(this, 0);
    public BigReal ToBigReal(int O) => BigReal.FromDouble(K, O);

    public Dble KMul(Rational k) => new(K * k);

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{K}";

    public static Dble NthRoot(Dble a, int n) => new(double.RootN(a.K, n));
    public static Dble Sqrt(Dble a) => new(double.Sqrt(a.K));

    public static Dble operator /(int a, Dble b) => a * b.Inv();

    public static double Abs(Dble t) => double.Abs(t.K);
    public static Dble Round(Dble e, int d, MidpointRounding mode = MidpointRounding.ToEven) 
        => new(double.Round(e.K, d, mode));

    public static bool IsValuedField => true;

    public static implicit operator double(Dble e) => e.K;

    public static implicit operator Dble(Rational e) => new(e);

    public static Dble operator +(Dble a, Rational b) => a.Add(b);

    public static Dble operator +(Rational a, Dble b) => b.Add(a);

    public static Dble operator -(Dble a, Rational b) => a.Sub(b);

    public static Dble operator -(Rational a, Dble b) => new Dble(a).Sub(b);

    public static Dble operator *(Dble a, Rational b) => a.Mul(new Dble(b));

    public static Dble operator *(Rational a, Dble b) => new Dble(a).Mul(b);

    public static Dble operator /(Dble a, Rational b) => a.Div(b).quo;

    public static bool operator ==(Dble a, Dble b) => Math.Abs(a.K - b.K) < EpsDouble;

    public static bool operator !=(Dble a, Dble b) => !(a == b);
    public static bool operator <(Dble a, Dble b) => a.K < b.K;

    public static bool operator >(Dble a, Dble b) => a.K > b.K;
    public static bool operator <=(Dble a, Dble b) => (a < b) || (a == b);

    public static bool operator >=(Dble a, Dble b) => (a > b) || (a == b);
    public static Dble Min(Dble a, Dble b) => a <= b ? a : b;
    public static Dble Max(Dble a, Dble b) => a < b ? b : a;
}