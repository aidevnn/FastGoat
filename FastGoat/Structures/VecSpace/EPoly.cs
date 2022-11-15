namespace FastGoat.Structures.VecSpace;

public readonly struct EPoly<K> : IVsElt<K, EPoly<K>>, IElt<EPoly<K>>, IRingElt<EPoly<K>>, IFieldElt<EPoly<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public KPoly<K> F { get; }
    public KPoly<K> Poly { get; }

    public EPoly(KPoly<K> f)
    {
        // if (f.P == 0)
        //     throw new GroupException(GroupExceptionType.GroupDef);

        F = f;
        Poly = f.X;
        Hash = (Poly.Hash, f.Hash).GetHashCode();
    }

    private EPoly(KPoly<K> f, KPoly<K> poly)
    {
        F = f;
        Poly = poly;
        Hash = (Poly.Hash, f.Hash).GetHashCode();
    }

    private EPoly(KPoly<K> f, K k)
    {
        F = f;
        Poly = new(F.x, k);
        Hash = (Poly.Hash, f.Hash).GetHashCode();
    }

    public bool Equals(EPoly<K> other) => Hash == other.Hash;

    public int CompareTo(EPoly<K> other) => Poly.CompareTo(other.Poly);

    public int P => F.P;
    public EPoly<K> Inv()
    {
        var (x, y) = Ring.Bezout(Poly, F);
        var gcd = (Poly * x + F * y).Div(F).rem;
        return new(F, (x / gcd).Div(F).rem);
    }

    public EPoly<K> KMul(K k) => new(F, Poly.KMul(k));

    public int Hash { get; }
    public bool IsZero() => Poly.IsZero();

    public K KZero => F.KZero;
    public K KOne => F.KOne;
    public EPoly<K> Zero => new(F, F.Zero);
    public EPoly<K> One => new(F, F.One);
    public EPoly<K> Add(EPoly<K> e)
    {
        return new(F, Poly.Add(e.Poly).Div(F).rem);
    }

    public EPoly<K> Sub(EPoly<K> e)
    {
        return new(F, Poly.Sub(e.Poly).Div(F).rem);
    }

    public EPoly<K> Opp() => new(F, Poly.Opp());

    public EPoly<K> Mul(EPoly<K> e) => new(F, Poly.Mul(e.Poly).Div(F).rem);

    public (EPoly<K> quo, EPoly<K> rem) Div(EPoly<K> e)
    {
        var (q, r) = Poly.Div(e.Poly);
        return (new(F, q.Div(F).rem), new(F, r.Div(F).rem));
    }

    public EPoly<K> Mul(int k) => new(F, Poly.Mul(k));

    public EPoly<K> Pow(int k) => new(F, Poly.Pow(k).Div(F).rem);

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";
        
        return $"{Poly}";
    }

    public static implicit operator EPoly<K>(KPoly<K> e) => new(e.One, e);
    
    public static EPoly<K> operator +(EPoly<K> a, EPoly<K> b) => a.Add(b);
    public static EPoly<K> operator -(EPoly<K> a, EPoly<K> b) => a.Sub(b);
    public static EPoly<K> operator -(EPoly<K> a) => a.Opp();
    public static EPoly<K> operator *(EPoly<K> a, EPoly<K> b) => a.Mul(b);
    public static EPoly<K> operator /(EPoly<K> a, EPoly<K> b) => a.Div(b).quo;

    public static EPoly<K> operator +(EPoly<K> a, K b) => a.Add(new(a.F, b));
    public static EPoly<K> operator +(K b, EPoly<K> a) => a.Add(new(a.F, b));
    public static EPoly<K> operator +(EPoly<K> a, int b) => a.Add(new(a.F, a.KOne.Mul(b)));
    public static EPoly<K> operator +(int b, EPoly<K> a) => a.Add(new(a.F, a.KOne.Mul(b)));

    public static EPoly<K> operator -(EPoly<K> a, K b) => a.Sub(new(a.F, b));
    public static EPoly<K> operator -(K b, EPoly<K> a) => new EPoly<K>(a.F, b).Sub(a);
    public static EPoly<K> operator -(EPoly<K> a, int b) => a.Sub(new(a.F, a.KOne.Mul(b)));
    public static EPoly<K> operator -(int b, EPoly<K> a) => new EPoly<K>(a.F, a.KOne.Mul(b)).Sub(a);

    public static EPoly<K> operator *(EPoly<K> a, K b) => a.KMul(b);
    public static EPoly<K> operator *(K b, EPoly<K> a) => a.KMul(b);
    public static EPoly<K> operator *(EPoly<K> a, int b) => a.Mul(b);
    public static EPoly<K> operator *(int b, EPoly<K> a) => a.Mul(b);

    public static EPoly<K> operator /(EPoly<K> a, K b) => a.KMul(b.Inv());
    public static EPoly<K> operator /(EPoly<K> a, int b) => a.KMul(a.KOne.Mul(b).Inv());
}