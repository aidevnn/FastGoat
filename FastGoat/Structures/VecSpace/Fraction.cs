namespace FastGoat.Structures.VecSpace;

public readonly struct Fraction<T> : IElt<Fraction<T>>, IRingElt<Fraction<T>>, IFieldElt<Fraction<T>>
    where T : struct, IElt<T>, IRingElt<T>
{
    public T Num { get; }
    public T Denom { get; }

    public Fraction(int p, T scalar)
    {
        P = p;

        // Console.WriteLine($"Check {typeof(T).GetProperties().Any(p => p.Name == "P")}"); // Reflection fashion

        Num = scalar;
        Denom = Num.One;
        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public Fraction(int p, T num, T denom)
    {
        P = p;
        if (denom.IsZero())
            throw new DivideByZeroException();

        var (q, r) = num.Div(denom);
        if (r.IsZero())
        {
            Num = q;
            Denom = q.One;
        }
        else
        {
            var gcd = Ring.Gcd(num, denom);
            var num0 = num.Div(gcd).quo;
            var denom0 = denom.Div(gcd).quo;
            var c = denom0.LeadingCoeff;
            Num = num0 / c;
            Denom = denom0 / c;
        }

        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public Fraction(int p, (T num, T denom) e)
    {
        P = p;
        Num = e.num;
        Denom = e.denom;
        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public Fraction<T> Zero => new(P, Num.Zero);
    public Fraction<T> One => new(P, Num.One);
    public Fraction<T> LeadingCoeff => One;
    public bool Equals(Fraction<T> other) => P == other.P && Num.Equals(other.Num) && Denom.Equals(other.Denom);

    public int CompareTo(Fraction<T> other) => (Num.Mul(other.Denom)).CompareTo(other.Num.Mul(Denom));

    public int P { get; }
    public Fraction<T> Inv() => new(P, Denom, Num);

    public int Hash { get; }
    public bool IsZero() => Num.IsZero();

    public Fraction<T> Add(Fraction<T> e) => new(P, Num.Mul(e.Denom).Add(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

    public Fraction<T> Sub(Fraction<T> e) => new(P, Num.Mul(e.Denom).Sub(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

    public Fraction<T> Opp() => new(P, Num.Opp(), Denom);

    public Fraction<T> Mul(Fraction<T> e) => new(P, Num.Mul(e.Num), Denom.Mul(e.Denom));

    public (Fraction<T> quo, Fraction<T> rem) Div(Fraction<T> e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        return (new(P, Num.Mul(e.Denom), Denom.Mul(e.Num)), Zero);
    }

    public Fraction<T> Mul(int k) => new(P, Num.Mul(k), Denom);

    public Fraction<T> Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var num = Num.Pow(k);
        var denom = Denom.Pow(k);
        return new(P, num, denom);
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";

        if (Denom.Equals(Denom.One))
            return $"{Num}";

        var num = $"{Num}".Contains('+') ? $"({Num})" : $"{Num}";
        var denom = $"{Denom}".Contains('+') ? $"({Denom})" : $"{Denom}";

        return $"{num} / {denom}";
    }

    public static Fraction<T> operator +(Fraction<T> a, Fraction<T> b) => a.Add(b);
    public static Fraction<T> operator +(int a, Fraction<T> b) => b.Add(b.One.Mul(a));
    public static Fraction<T> operator +(Fraction<T> a, int b) => a.Add(a.One.Mul(b));
    public static Fraction<T> operator -(Fraction<T> a) => a.Opp();
    public static Fraction<T> operator -(Fraction<T> a, Fraction<T> b) => a + (-b);
    public static Fraction<T> operator -(int a, Fraction<T> b) => a + (-b);
    public static Fraction<T> operator -(Fraction<T> a, int b) => a + (-b);
    public static Fraction<T> operator *(Fraction<T> a, Fraction<T> b) => a.Mul(b);
    public static Fraction<T> operator *(Fraction<T> a, int b) => a.Mul(b);
    public static Fraction<T> operator *(int a, Fraction<T> b) => b.Mul(a);
    public static Fraction<T> operator /(Fraction<T> a, Fraction<T> b) => a.Div(b).quo;
    public static Fraction<T> operator /(Fraction<T> a, int b) => a.Div(a.One.Mul(b)).quo;
    public static Fraction<T> operator /(int a, Fraction<T> b) => b.Inv().Mul(a);
}