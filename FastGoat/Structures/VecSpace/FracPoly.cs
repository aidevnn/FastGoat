namespace FastGoat.Structures.VecSpace;

public readonly struct FracPoly<K> : IVsElt<K, FracPoly<K>>, IElt<FracPoly<K>>, IRingElt<FracPoly<K>>,
    IFieldElt<FracPoly<K>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public KPoly<K> Num { get; }
    public KPoly<K> Denom { get; }

    public FracPoly(char x, K scalar)
    {
        var num = new KPoly<K>(x, scalar);
        Num = num.X;
        P = Num.P;
        Denom = Num.One;
        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public static bool IsValuedField => false;
    public static double Abs(FracPoly<K> e) => throw new NotImplementedException();

    public FracPoly(KPoly<K> num)
    {
        Num = num;
        P = Num.P;
        Denom = Num.One;
        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public FracPoly(KPoly<K> num, KPoly<K> denom)
    {
        if (num.P != denom.P || num.x != denom.x)
            throw new GroupException(GroupExceptionType.GroupDef);

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
            var gcd = Ring.FastGCD(num, denom);
            var num0 = num.Div(gcd).quo;
            var denom0 = denom.Div(gcd).quo;
            var c = denom0.Coefs.Last();
            Num = num0 / c;
            Denom = denom0 / c;
        }

        Hash = (Num.Hash, Denom.Hash).GetHashCode();
    }

    public char x => Num.x;

    public bool Equals(FracPoly<K> other) => Num.Equals(other.Num) && Denom.Equals(other.Denom);

    public int CompareTo(FracPoly<K> other) => (Num * other.Denom).CompareTo(Denom * other.Num);

    public K KZero => Num.KZero;
    public K KOne => Num.KOne;

    public FracPoly<K> KMul(K k) => new(k * Num, Denom);
    public int Hash { get; }

    public bool IsZero() => Num.IsZero();

    public FracPoly<K> Substitute(FracPoly<K> s)
    {
        var num = Num.Coefs.Select((c, i) => c * s.Pow(i)).Aggregate(Zero, (acc, a) => a + acc);
        var denom = Denom.Coefs.Select((c, i) => c * s.Pow(i)).Aggregate(Zero, (acc, a) => a + acc);
        return num / denom;
    }

    public FracPoly<K> Substitute(KPoly<K> s)
    {
        var num = Num.Coefs.Select((c, i) => c * s.Pow(i)).Aggregate(s.Zero, (acc, a) => a + acc);
        var denom = Denom.Coefs.Select((c, i) => c * s.Pow(i)).Aggregate(s.Zero, (acc, a) => a + acc);
        return new(num, denom);
    }

    public int P { get; }
    public FracPoly<K> Zero => new(Num.Zero);
    public FracPoly<K> One => new(Num.One);
    public FracPoly<K> X => new(Num.X);

    public FracPoly<K> Add(FracPoly<K> e) => new(Num.Mul(e.Denom).Add(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

    public FracPoly<K> Sub(FracPoly<K> e) => new(Num.Mul(e.Denom).Sub(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

    public FracPoly<K> Opp() => new(Num.Opp(), Denom);

    public FracPoly<K> Mul(FracPoly<K> e) => new(Num * e.Num, Denom * e.Denom);

    public FracPoly<K> Inv() => new(Denom, Num);
    public bool Invertible() => !Num.IsZero();

    public (FracPoly<K> quo, FracPoly<K> rem) Div(FracPoly<K> e) =>
        (new FracPoly<K>(Num * e.Denom, Denom * e.Num), Zero);

    public FracPoly<K> Mul(int k) => new(k * Num, Denom);

    public FracPoly<K> Pow(int k)
    {
        if (k < 0)
            return new(Denom.Pow(-k), Num.Pow(-k));
        
        return new(Num.Pow(k), Denom.Pow(k));
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (Denom.Equals(Denom.One))
            return Num.ToString();

        return $"({Num})/({Denom})";
    }

    public static FracPoly<K> operator +(FracPoly<K> a, FracPoly<K> b) => a.Add(b);

    public static FracPoly<K> operator +(int a, FracPoly<K> b) => (b.One.Mul(a)).Add(b);

    public static FracPoly<K> operator +(FracPoly<K> a, int b) => a.Add(a.One.Mul(b));

    public static FracPoly<K> operator -(FracPoly<K> a) => a.Opp();

    public static FracPoly<K> operator -(FracPoly<K> a, FracPoly<K> b) => a.Sub(b);

    public static FracPoly<K> operator -(int a, FracPoly<K> b) => (b.One.Mul(a)).Sub(b);

    public static FracPoly<K> operator -(FracPoly<K> a, int b) => a.Sub(a.One.Mul(b));

    public static FracPoly<K> operator *(FracPoly<K> a, FracPoly<K> b) => a.Mul(b);

    public static FracPoly<K> operator *(int a, FracPoly<K> b) => b.Mul(a);

    public static FracPoly<K> operator *(FracPoly<K> a, int b) => a.Mul(b);

    public static FracPoly<K> operator /(FracPoly<K> a, FracPoly<K> b) => a.Div(b).quo;

    public static FracPoly<K> operator /(FracPoly<K> a, int b) => new(a.Num, a.Denom * b);

    public static FracPoly<K> operator /(int a, FracPoly<K> b) => new(b.Denom, a * b.Num);

    public static FracPoly<K> operator +(FracPoly<K> a, K b) => a + new FracPoly<K>(new KPoly<K>(a.x, b));

    public static FracPoly<K> operator +(K a, FracPoly<K> b) => new FracPoly<K>(new KPoly<K>(b.x, a)) + b;

    public static FracPoly<K> operator -(FracPoly<K> a, K b) => a - new FracPoly<K>(new KPoly<K>(a.x, b));

    public static FracPoly<K> operator -(K a, FracPoly<K> b) => new FracPoly<K>(new KPoly<K>(b.x, a)) - b;

    public static FracPoly<K> operator *(FracPoly<K> a, K b) => new(a.Num * b, a.Denom);

    public static FracPoly<K> operator *(K a, FracPoly<K> b) => new(a * b.Num, b.Denom);

    public static FracPoly<K> operator /(FracPoly<K> a, K b) => new(a.Num, a.Denom * b);
}