namespace FastGoat.Structures.VecSpace;

public readonly struct EPolynomial<K> : IVsElt<K, EPolynomial<K>>, IElt<EPolynomial<K>>,
    IRingElt<EPolynomial<K>>, IFieldElt<EPolynomial<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public Polynomial<K, Xi> Num { get; }
    public Polynomial<K, Xi> Denom { get; }
    public PolynomialBasis<K, Xi> Basis { get; }
    public int P => Num.P;

    public static bool IsValuedField => false;
    public static double Abs(EPolynomial<K> e) => throw new NotImplementedException();

    public EPolynomial(Polynomial<K, Xi> num, PolynomialBasis<K, Xi> basis)
    {
        if (!num.Indeterminates.Equals(basis.Indeterminates))
            throw new GroupException(GroupExceptionType.GroupDef);

        Num = num;
        Denom = num.One;
        Basis = basis;
        Hash = Basis.Hash;
    }

    public EPolynomial(Polynomial<K, Xi> num, Polynomial<K, Xi> denom)
    {
        if (num.P != denom.P || !num.Indeterminates.Equals(denom.Indeterminates))
            throw new GroupException(GroupExceptionType.GroupDef);

        if (denom.IsZero())
            throw new DivideByZeroException();

        Basis = new PolynomialBasis<K, Xi>(num.Indeterminates);
        (Num, Denom) = SimplifyFrac(Basis, num, denom);
        Hash = Basis.Hash;
    }

    public EPolynomial(Polynomial<K, Xi> num, Polynomial<K, Xi> denom, PolynomialBasis<K, Xi> basis)
    {
        Basis = basis;
        (Num, Denom) = SimplifyFrac(Basis, num, denom);
        Hash = Basis.Hash;
    }

    public int Hash { get; }
    public bool Equals(EPolynomial<K> other) => (Num * other.Denom).Equals(Denom * other.Num);

    public int CompareTo(EPolynomial<K> other) => (Num * other.Denom).CompareTo(Denom * other.Num);
    public K KZero => Num.KZero;
    public K KOne => Denom.KOne;
    public bool IsZero() => Num.IsZero();

    public EPolynomial<K> Zero => new(Num.Zero, Basis);
    public EPolynomial<K> One => new(Num.One, Basis);

    public EPolynomial<K> Add(EPolynomial<K> e)
    {
        var num = Basis.Rem(Num * e.Denom + Denom * e.Num);
        var denom = Basis.Rem(Denom * e.Denom);
        return new(num, denom, Basis);
    }

    public EPolynomial<K> Sub(EPolynomial<K> e)
    {
        var num = Basis.Rem(Num * e.Denom - Denom * e.Num);
        var denom = Basis.Rem(Denom * e.Denom);
        return new(num, denom, Basis);
    }

    public EPolynomial<K> Opp()
    {
        var num = Basis.Rem(-Num);
        return new(num, Denom, Basis);
    }

    public EPolynomial<K> Mul(EPolynomial<K> e)
    {
        var num = Basis.Rem(Num * e.Num);
        var denom = Basis.Rem(Denom * e.Denom);
        return new(num, denom, Basis);
    }

    public (EPolynomial<K> quo, EPolynomial<K> rem) Div(EPolynomial<K> e)
    {
        var num = Basis.Rem(Num * e.Denom);
        var denom = Basis.Rem(Denom * e.Num);
        return (new(num, denom, Basis), Zero);
    }

    public EPolynomial<K> Inv() => new(Denom, Num, Basis);
    public bool Invertible() => !IsZero();

    public EPolynomial<K> Mul(int k)
    {
        var num = Basis.Rem(k * Num);
        return new(num, Denom, Basis);
    }

    public EPolynomial<K> Pow(int k)
    {
        var num = Basis.Rem(Num.Pow(k));
        var denom = Basis.Rem(Denom.Pow(k));
        return new(num, denom, Basis);
    }

    public EPolynomial<K> KMul(K k)
    {
        var num = Basis.Rem(k * Num);
        return new(num, Denom, Basis);
    }

    public override string ToString()
    {
        var num = Num.ToString();
        var denom = Denom.Equals(Denom.One) ? "" : Denom.ToString();
        num = num.Contains('+') ? $"({num})" : num;
        denom = denom.Contains('+') ? $"({denom})" : denom;

        if (denom == "")
            return num;

        return $"{num}/{denom}";
    }

    public override int GetHashCode() => Hash;

    public static EPolynomial<K> operator +(EPolynomial<K> a, EPolynomial<K> b) => a.Add(b);
    public static EPolynomial<K> operator +(int a, EPolynomial<K> b) => b.Add(b.One.Mul(a));
    public static EPolynomial<K> operator +(EPolynomial<K> a, int b) => a.Add(a.One.Mul(b));
    public static EPolynomial<K> operator -(EPolynomial<K> a) => a.Opp();
    public static EPolynomial<K> operator -(EPolynomial<K> a, EPolynomial<K> b) => a + (-b);
    public static EPolynomial<K> operator -(int a, EPolynomial<K> b) => a + (-b);
    public static EPolynomial<K> operator -(EPolynomial<K> a, int b) => a + (-b);
    public static EPolynomial<K> operator *(EPolynomial<K> a, EPolynomial<K> b) => a.Mul(b);
    public static EPolynomial<K> operator *(EPolynomial<K> a, int b) => a.Mul(b);
    public static EPolynomial<K> operator *(int a, EPolynomial<K> b) => b.Mul(a);
    public static EPolynomial<K> operator /(EPolynomial<K> a, EPolynomial<K> b) => a.Div(b).quo;
    public static EPolynomial<K> operator /(EPolynomial<K> a, int b) => a.Div(a.One.Mul(b)).quo;
    public static EPolynomial<K> operator /(int a, EPolynomial<K> b) => b.Inv() * a;
    public static EPolynomial<K> operator +(EPolynomial<K> a, K b) => a + a.One.KMul(b);
    public static EPolynomial<K> operator +(K a, EPolynomial<K> b) => b.One.KMul(a) + b;
    public static EPolynomial<K> operator -(EPolynomial<K> a, K b) => a - a.One.KMul(b);
    public static EPolynomial<K> operator -(K a, EPolynomial<K> b) => b.One.KMul(a) - b;
    public static EPolynomial<K> operator *(EPolynomial<K> a, K b) => a.KMul(b);
    public static EPolynomial<K> operator *(K a, EPolynomial<K> b) => b.KMul(a);
    public static EPolynomial<K> operator /(EPolynomial<K> a, K b) => a.KMul(b.Inv());

    public static (Polynomial<K, Xi> Num, Polynomial<K, Xi> Denom) SimplifyFrac(PolynomialBasis<K, Xi> basis, Polynomial<K, Xi> num,
        Polynomial<K, Xi> denom)
    {
        var (q, r) = num.Div(denom);
        var r0 = basis.Rem(r);
        if (r0.IsZero())
        {
            var q0 = basis.Rem(q);
            return (q0, q0.One);
        }

        var nbNum = num.NbIndeterminates;
        var nbDenom = denom.NbIndeterminates;
        if (nbNum == 0)
            return (num, denom);

        if (nbNum == 1 && nbDenom == 1)
        {
            var ai = num.ExtractIndeterminate;
            var bi = denom.ExtractIndeterminate;
            if (!ai.Equals(bi))
                return (num, denom);

            var num1 = num.ToKPoly(ai);
            var denom1 = denom.ToKPoly(bi);
            var gcd = Ring.Gcd(num1, denom1).Monic;
            var num2 = basis.Rem((num1 / gcd).ToPolynomial(num.Indeterminates, ai));
            var denom2 = basis.Rem((denom1 / gcd).ToPolynomial(num.Indeterminates, ai));
            return (num2, denom2);
        }

        if (nbNum > 3 || nbDenom > 3)
            return (num, denom);

        var lcm = Ring.LcmPolynomial(num, denom);
        var num0 = basis.Rem(lcm / denom);
        var denom0 = basis.Rem(lcm / num);
        return (num0, denom0);
    }
}