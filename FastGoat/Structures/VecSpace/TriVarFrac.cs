namespace FastGoat.Structures.VecSpace;

public readonly struct TriVarFrac<K> : IElt<TriVarFrac<K>>, IRingElt<TriVarFrac<K>>, IFieldElt<TriVarFrac<K>>,
    IModuleElt<K, TriVarFrac<K>>, IVsElt<K, TriVarFrac<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public TriVarPoly<K> Num { get; }
    public TriVarPoly<K> Denom { get; }
    public TriVarFracSimplifier<K> Simplifier { get; }

    public TriVarFrac(TriVarFracSimplifier<K> simp, TriVarPoly<K> num, TriVarPoly<K> denom)
    {
        (Simplifier, Num, Denom) = (simp, num, denom);

        if (denom.IsZero())
            throw new DivideByZeroException();
        else
        {
            var (q, r) = num.Div(denom);
            if (r.IsZero())
                (Num, Denom) = (q, q.One);
        }

        var lt = Denom.LeadingDetails.lc;
        if (!lt.IsZero())
        {
            Num *= lt.Inv();
            Denom *= lt.Inv();
        }

        Hash = (Num.Hash, Denom.Hash, Simplifier.GetHashCode()).GetHashCode();
    }

    public TriVarFrac(TriVarPoly<K> num, TriVarPoly<K> denom) : this(new TriVarFracSimplifier<K>(), num, denom)
    {
    }

    public bool Equals(TriVarFrac<K> other)
    {
        return Sub(other).IsZero();
    }

    public int CompareTo(TriVarFrac<K> other)
    {
        return (Num * other.Denom).CompareTo(Denom * other.Num);
    }

    public int Hash { get; }
    public bool IsZero() => !Denom.IsZero() && Num.IsZero();

    public bool IsDivZero() => Simplifier.IsDivZero(this);

    public TriVarFrac<K> Zero => new(Simplifier, Num.Zero, Num.One);
    public TriVarFrac<K> One => new(Simplifier, Num.One, Num.One);

    public K KZero => Num.KZero;
    public K KOne => Num.KOne;
    public TriVarFrac<K> X1 => new(Simplifier, Num.X1, Num.One);
    public TriVarFrac<K> X2 => new(Simplifier, Num.X2, Num.One);
    public TriVarFrac<K> X3 => new(Simplifier, Num.X3, Num.One);

    public TriVarFrac<K> KMul(K k)
    {
        return new(Simplifier, k * Num, Denom);
    }

    public TriVarFrac<K> Add(TriVarFrac<K> e)
    {
        var num = (Num * e.Denom + Denom * e.Num);
        var denom = Denom.Mul(e.Denom);
        return Simplifier.Apply(num, denom);
    }

    public TriVarFrac<K> Sub(TriVarFrac<K> e)
    {
        var num = (Num * e.Denom - Denom * e.Num);
        var denom = Denom.Mul(e.Denom);
        return Simplifier.Apply(num, denom);
    }

    public TriVarFrac<K> Opp() => new(Simplifier, -Num, Denom);

    public TriVarFrac<K> Mul(TriVarFrac<K> e)
    {
        var num = (Num * e.Num);
        var denom = (Denom * e.Denom);
        return Simplifier.Apply(num, denom);
    }

    public (TriVarFrac<K> quo, TriVarFrac<K> rem) Div(TriVarFrac<K> e)
    {
        var num = Num * e.Denom;
        var denom = Denom * e.Num;
        return (Simplifier.Apply(num, denom), Zero);
    }

    public TriVarFrac<K> Mul(int k) => KMul(k * KOne);

    public TriVarFrac<K> Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        return this.FastPow(k);
    }

    public int P => Num.P;

    public TriVarFrac<K> Inv()
    {
        return new(Simplifier, Denom, Num);
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (Denom.Equals(Denom.One))
            return Num.ToString();

        return $"({Num})/({Denom})";
    }

    public bool Invertible() => !IsDivZero();

    public static TriVarFrac<K> operator +(TriVarFrac<K> a, TriVarFrac<K> b) => a.Add(b);

    public static TriVarFrac<K> operator +(int a, TriVarFrac<K> b) => (b.One.Mul(a)).Add(b);

    public static TriVarFrac<K> operator +(TriVarFrac<K> a, int b) => a.Add(a.One.Mul(b));

    public static TriVarFrac<K> operator -(TriVarFrac<K> a) => a.Opp();

    public static TriVarFrac<K> operator -(TriVarFrac<K> a, TriVarFrac<K> b) => a.Sub(b);

    public static TriVarFrac<K> operator -(int a, TriVarFrac<K> b) => (b.One.Mul(a)).Sub(b);

    public static TriVarFrac<K> operator -(TriVarFrac<K> a, int b) => a.Sub(a.One.Mul(b));

    public static TriVarFrac<K> operator *(TriVarFrac<K> a, TriVarFrac<K> b) => a.Mul(b);

    public static TriVarFrac<K> operator *(int a, TriVarFrac<K> b) => b.Mul(a);

    public static TriVarFrac<K> operator *(TriVarFrac<K> a, int b) => a.Mul(b);

    public static TriVarFrac<K> operator /(TriVarFrac<K> a, TriVarFrac<K> b) => a.Div(b).quo;

    public static TriVarFrac<K> operator /(TriVarFrac<K> a, int b) => new(a.Simplifier, a.Num, a.Denom * b);

    public static TriVarFrac<K> operator /(int a, TriVarFrac<K> b) => a * b.Inv();

    public static TriVarFrac<K> operator +(TriVarFrac<K> a, K b) => a + a.One.KMul(b);

    public static TriVarFrac<K> operator +(K a, TriVarFrac<K> b) => b.One.KMul(a) + b;

    public static TriVarFrac<K> operator -(TriVarFrac<K> a, K b) => a - a.One.KMul(b);

    public static TriVarFrac<K> operator -(K a, TriVarFrac<K> b) => b.One.KMul(a) - b;

    public static TriVarFrac<K> operator *(TriVarFrac<K> a, K b) => a.KMul(b);

    public static TriVarFrac<K> operator *(K a, TriVarFrac<K> b) => b.KMul(a);

    public static TriVarFrac<K> operator /(TriVarFrac<K> a, K b) => new(a.Simplifier, a.Num, a.Denom * b);
    public static double Abs(TriVarFrac<K> t) => TriVarPoly<K>.Abs(t.Num) / TriVarPoly<K>.Abs(t.Denom);

    public static bool IsValuedField => false;
}