using System.Collections;

namespace FastGoat.Structures.VecSpace;

public readonly struct EPoly<K> : IElt<EPoly<K>>, IRingElt<EPoly<K>>, IFieldElt<EPoly<K>>, IModuleElt<K, EPoly<K>>,
    IVsElt<K, EPoly<K>>, IEnumerable<K>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public KPoly<K> F { get; }
    public KPoly<K> Poly { get; }
    public static double Abs(EPoly<K> e) => throw new NotImplementedException();
    public static bool IsValuedField => false;

    public EPoly(KPoly<K> f)
    {
        // if (f.P == 0)
        //     throw new GroupException(GroupExceptionType.GroupDef);

        F = f;
        Poly = f.X;
        // Hash = (Poly.Hash, f.Hash).GetHashCode();
        Hash = (Poly.Degree, f.Degree).GetHashCode();
    }

    public EPoly(KPoly<K> f, KPoly<K> poly)
    {
        F = f;
        Poly = poly;
        // Hash = (Poly.Hash, f.Hash).GetHashCode();
        Hash = (Poly.Degree, f.Degree).GetHashCode();
    }

    public int GetHashCodeSlow() => (F.GetHashCodeSlow(), Poly.GetHashCodeSlow()).GetHashCode();
    public bool Equals(EPoly<K> other) => Poly.Equals(other.Poly); // Avoid collisions

    public int CompareTo(EPoly<K> other) => Poly.CompareTo(other.Poly);

    public int P => F.P;
    public EPoly<K> Clone => new(F, Poly);

    public EPoly<K> Inv()
    {
        var (x, y) = P == 0 ? Ring.FastBezout(Poly, F) : Ring.Bezout(Poly, F);
        var gcd = (Poly * x + F * y).Div(F).rem;
        return new(F, (x / gcd).Div(F).rem);
    }

    public IEnumerator<K> GetEnumerator() => Poly.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public bool Invertible() => !IsZero();

    public EPoly<K> KMul(K k) => new(F, Poly.KMul(k).Div(F).rem);

    public int Hash { get; }
    public bool IsZero() => Poly.IsZero();
    public K this[int idx] => Poly[idx];
    public int Degree => Poly.Degree;

    public K KZero => F.KZero;
    public K KOne => F.KOne;
    public EPoly<K> Zero => new(F, F.Zero);
    public EPoly<K> One => new(F, F.One);
    public EPoly<K> X => new(F, F.X.Div(F).rem);
    public EPoly<K> Derivative => new(F, Poly.Derivative.Div(F).rem);
    public EPoly<K> Substitute(KPoly<K> f) => new(F, Poly.Substitute(f).Div(F).rem);
    public EPoly<K> Substitute(EPoly<K> f) => Poly.Substitute(f);

    public FracPoly<K> Substitute(FracPoly<K> f)
    {
        return Poly.Coefs.Select((k, i) => k * f.Pow(i)).Aggregate(f.Zero, (acc, a) => acc + a);
    }

    public EPoly<K> Compose(int k)
    {
        if (k < 0)
            throw new();

        if (k == 0)
            return X;

        return Substitute(Compose(k - 1));
    }

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
        var gcd = P == 0 ? Ring.FastGCD(Poly, e.Poly).Monic : Ring.Gcd(Poly, e.Poly);
        var e0 = new EPoly<K>(F, Poly / gcd);
        var e1 = new EPoly<K>(F, e.Poly / gcd);
        return (e0.Mul(e1.Inv()), Zero);
    }

    public EPoly<K> Mul(int k) => new(F, Poly.Mul(k).Div(F).rem);

    public EPoly<K> Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        if (Degree == 0)
            return new(F, Poly.Pow(k));

        return Ring.FastPow(this, k);
    }

    public KPoly<EPoly<K>> ToKPoly(char x) => new(x, this);
    public FracPoly<EPoly<K>> ToFracPoly(char x) => new(ToKPoly(x));

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";

        return $"{Poly}";
    }

    public static EPoly<K> operator +(EPoly<K> a, EPoly<K> b) => a.Add(b);
    public static EPoly<K> operator +(int a, EPoly<K> b) => b.Add(b.One.Mul(a));
    public static EPoly<K> operator +(EPoly<K> a, int b) => a.Add(a.One.Mul(b));
    public static EPoly<K> operator -(EPoly<K> a) => a.Opp();
    public static EPoly<K> operator -(EPoly<K> a, EPoly<K> b) => a + (-b);
    public static EPoly<K> operator -(int a, EPoly<K> b) => a + (-b);
    public static EPoly<K> operator -(EPoly<K> a, int b) => a + (-b);
    public static EPoly<K> operator *(EPoly<K> a, EPoly<K> b) => a.Mul(b);
    public static EPoly<K> operator *(EPoly<K> a, int b) => a.Mul(b);
    public static EPoly<K> operator *(int a, EPoly<K> b) => b.Mul(a);
    public static EPoly<K> operator /(EPoly<K> a, EPoly<K> b) => a.Div(b).quo;
    public static EPoly<K> operator /(EPoly<K> a, int b) => a.Div(a.One.Mul(b)).quo;
    public static EPoly<K> operator /(int a, EPoly<K> b) => b.Inv().Mul(a);

    public static EPoly<K> operator +(EPoly<K> a, K b) => a + a.One.KMul(b);
    public static EPoly<K> operator +(K a, EPoly<K> b) => b.One.KMul(a) + b;
    public static EPoly<K> operator -(EPoly<K> a, K b) => a - a.One.KMul(b);
    public static EPoly<K> operator -(K a, EPoly<K> b) => b.One.KMul(a) - b;
    public static EPoly<K> operator *(EPoly<K> a, K b) => a.KMul(b);
    public static EPoly<K> operator *(K a, EPoly<K> b) => b.KMul(a);
    public static EPoly<K> operator /(EPoly<K> a, K b) => a.KMul(b.Inv());
    public static EPoly<K> operator /(K a, EPoly<K> b) => b.Inv().KMul(a);
    public static EPoly<K> operator +(EPoly<K> a, KPoly<K> b) => a + new EPoly<K>(a.F, b.Div(a.F).rem);
    public static EPoly<K> operator +(KPoly<K> a, EPoly<K> b) => new EPoly<K>(b.F, a.Div(b.F).rem) + b;
    public static EPoly<K> operator -(EPoly<K> a, KPoly<K> b) => a - new EPoly<K>(a.F, b.Div(a.F).rem);
    public static EPoly<K> operator -(KPoly<K> a, EPoly<K> b) => new EPoly<K>(b.F, a.Div(b.F).rem) - b;
    public static EPoly<K> operator *(EPoly<K> a, KPoly<K> b) => a * new EPoly<K>(a.F, b.Div(a.F).rem);
    public static EPoly<K> operator *(KPoly<K> a, EPoly<K> b) => new EPoly<K>(b.F, a.Div(b.F).rem) * b;
    public static EPoly<K> operator /(EPoly<K> a, KPoly<K> b) => a.Mul(new EPoly<K>(a.F, b.Div(a.F).rem).Inv());
    public static EPoly<K> operator /(KPoly<K> a, EPoly<K> b) => new EPoly<K>(b.F, a.Div(b.F).rem) * b.Inv();
}