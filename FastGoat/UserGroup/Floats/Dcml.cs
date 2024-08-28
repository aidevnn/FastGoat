using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.Floats;

public readonly struct Dcml : IElt<Dcml>, IRingElt<Dcml>, IFieldElt<Dcml>, IVsElt<Rational, Dcml>,
    IFloatElt<Dcml>, IFixedPrecisionElt<Dcml>
{
    public override bool Equals(object? obj)
    {
        return obj is Dcml other && Equals(other);
    }

    public static double Eps => 1.0e-26;
    public static Dcml Pi(int o = 50) => 3.14159265358979323846264338327m;

    public static Dcml DcmlZero() => new(0.0m);
    public static Dcml DcmlOne() => new(1.0m);
    public decimal K { get; }

    public Dcml(decimal k)
    {
        K = decimal.Abs(k) < (decimal)(Eps * 10) ? 0 : k;
        Hash = K.GetHashCode();
    }

    public bool Equals(Dcml other) => decimal.Abs(K - other.K) < (decimal)(Eps * 100);

    public int CompareTo(Dcml other) => K.CompareTo(other.K);

    public int Hash { get; }
    public bool IsZero() => decimal.Abs(K) < (decimal)(Eps * 100);

    public Dcml Zero => new(0);
    public Dcml One => new(1);
    public Dcml Add(Dcml e) => new(K + e.K);

    public Dcml Sub(Dcml e) => new(K - e.K);

    public Dcml Opp() => new(-K);

    public Dcml Mul(Dcml e) => new(K * e.K);

    public (Dcml quo, Dcml rem) Div(Dcml e) => (new(K / e.K), Zero);

    public Dcml Mul(int k) => new(K * k);

    public Dcml Pow(int k)
    {
        if (k < 0)
            return Inv().Pow(-k);

        var d0 = 1.0m;
        for (int i = 0; i < k; ++i)
            d0 *= K;

        return new(d0);
    }

    public static Dcml operator +(Dcml a, Dcml b) => a.Add(b);

    public static Dcml operator +(int a, Dcml b) => new(a + b.K);

    public static Dcml operator +(Dcml a, int b) => new(a.K + b);

    public static Dcml operator -(Dcml a) => a.Opp();

    public static Dcml operator -(Dcml a, Dcml b) => a.Sub(b);

    public static Dcml operator -(int a, Dcml b) => new(a - b.K);

    public static Dcml operator -(Dcml a, int b) => a.Sub(new(b));

    public static Dcml operator *(Dcml a, Dcml b) => a.Mul(b);

    public static Dcml operator *(int a, Dcml b) => b.Mul(a);

    public static Dcml operator *(Dcml a, int b) => a.Mul(b);

    public static Dcml operator /(Dcml a, Dcml b) => a.Div(b).quo;

    public static Dcml operator /(Dcml a, int b) => new(a.K / b);

    public int P => 0;
    public Dcml Inv() => new(1.0m / K);
    public int Sign => decimal.Sign(K);

    public bool Invertible() => !IsZero();
    public Rational KZero => Rational.KZero();
    public Rational KOne => Rational.KOne();
    public Dcml Absolute => new(Sign * K);
    public static Dcml From<T>(T e) where T : IElt<T>, IRingElt<T>, IFieldElt<T>, IFloatElt<T>
    {
        if (e is BigReal e0)
            return e0.ToDcml;

        if (e is Dble e1)
            return new((decimal)e1.K);

        if (e is Rational e2)
            return new((decimal)e2.ToDouble); // TODO missing digits

        if (e is Dcml e3)
            return e3;

        throw new ArgumentException();
    }

    public static int Digits => 26;

    public Dcml RoundEven => Round(this, 0);
    public BigReal ToBigReal(int O) => BigReal.FromString($"{K:E28}", O);

    public Dcml KMul(Rational k) => new(K * (decimal)k.ToDouble);

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{K}";

    public static Dcml NthRoot(Dcml r, int n)
    {
        if (n == 0)
            return r.One;

        if (n < 0)
            return NthRoot(r, -n).Inv();

        if (n % 2 == 0 && r.Sign == -1)
            throw new("Even NthRoot must has positive argument");

        if (n == 1)
            return r;

        Dcml ai;
        var aj = 1 + r / n;
        do
        {
            ai = aj;
            var aiPow = ai.Pow(n - 1);
            var num = aiPow * ai - r;
            var denom = n * aiPow;
            aj = ai - num / denom; // Newton iteration
        } while (!(ai - aj).IsZero());

        return aj;
    }
    public static Dcml Sqrt(Dcml a) => NthRoot(a, 2);
    public double ToDouble => (double)K;

    public static Dcml operator /(int a, Dcml b) => a * b.Inv();

    public static double Abs(Dcml t) => (double)decimal.Abs(t.K);
    public static Dcml Round(Dcml e, int d, MidpointRounding mode = MidpointRounding.ToEven) 
        => new(decimal.Round(e.K, d, mode));

    public static bool IsValuedField => true;

    public static implicit operator decimal(Dcml e) => e.K;
    public static implicit operator Dcml(decimal e) => new(e);

    public static implicit operator Dcml(Rational e) => new((decimal)e.ToDouble);

    public static Dcml operator +(Dcml a, Rational b) => a.Add(b);

    public static Dcml operator +(Rational a, Dcml b) => b.Add(a);

    public static Dcml operator -(Dcml a, Rational b) => a.Sub(b);

    public static Dcml operator -(Rational a, Dcml b) => new Dcml((decimal)a.ToDouble).Sub(b);

    public static Dcml operator *(Dcml a, Rational b) => a.Mul(new Dcml((decimal)b.ToDouble));

    public static Dcml operator *(Rational a, Dcml b) => new Dcml((decimal)a.ToDouble).Mul(b);

    public static Dcml operator /(Dcml a, Rational b) => a.Div(b).quo;

    public static bool operator ==(Dcml a, Dcml b) => a.Equals(b);

    public static bool operator !=(Dcml a, Dcml b) => !(a == b);
    public static bool operator <(Dcml a, Dcml b) => a.K < b.K;

    public static bool operator >(Dcml a, Dcml b) => a.K > b.K;
    public static bool operator <=(Dcml a, Dcml b) => (a < b) || (a == b);

    public static bool operator >=(Dcml a, Dcml b) => (a > b) || (a == b);
    public static Dcml Min(Dcml a, Dcml b) => a <= b ? a : b;
    public static Dcml Max(Dcml a, Dcml b) => a < b ? b : a;
}