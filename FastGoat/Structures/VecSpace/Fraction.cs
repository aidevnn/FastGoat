namespace FastGoat.Structures.VecSpace;

public readonly struct Fraction<T> : IElt<Fraction<T>>, IRingElt<Fraction<T>>, IFieldElt<Fraction<T>>
    where T : struct, IElt<T>, IRingElt<T>
{
    public T Num { get; }
    public T Denom { get; }

    public Fraction(T scalar)
    {
        P = 0;
        Num = scalar;
        Denom = Num.One;
        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public Fraction(T num, T denom)
    {
        P = 0;
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
            Num = num.Div(gcd).quo;
            Denom = denom.Div(gcd).quo;
        }
        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public Fraction<T> Zero => new(Num.Zero);
    public Fraction<T> One => new(Num.One);

    public bool Equals(Fraction<T> other) => Num.Equals(other.Num) && Denom.Equals(other.Denom);

    public int CompareTo(Fraction<T> other) => (Num.Mul(other.Denom)).CompareTo(other.Num.Mul(Denom));

    public int P { get; }
    public Fraction<T> Inv() => new(Denom, Num);

    public int Hash { get; }
    public bool IsZero() => Num.IsZero();

    public Fraction<T> Add(Fraction<T> e) => new(Num.Mul(e.Denom).Add(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

    public Fraction<T> Sub(Fraction<T> e) => new(Num.Mul(e.Denom).Sub(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

    public Fraction<T> Opp() => new(Num.Opp(), Denom);

    public Fraction<T> Mul(Fraction<T> e) => new(Num.Mul(e.Num), Denom.Mul(e.Denom));

    public (Fraction<T> quo, Fraction<T> rem) Div(Fraction<T> e)
    {
        if (e.IsZero())
            throw new DivideByZeroException();

        return (new(Num.Mul(e.Denom), Denom.Mul(e.Num)), Zero);
    }

    public Fraction<T> Mul(int k) => new(Num.Mul(k), Denom);

    public Fraction<T> Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        var num = Num.Pow(k);
        var denom = Denom.Pow(k);
        return new(num, denom);
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