using System.Collections;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;

namespace FastGoat.Structures.VecSpace;

public readonly struct SPoly<K> : IElt<SPoly<K>>, IRingElt<SPoly<K>>, IFieldElt<SPoly<K>>, IModuleElt<K, SPoly<K>>,
    IVsElt<K, SPoly<K>>, IEnumerable<K>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public int Ord { get; }
    public KPoly<K> Poly { get; }

    public static double Abs(SPoly<K> e) => e.IsZero()
        ? double.PositiveInfinity
        : e.Index().First(c => !c.Item.IsZero()).Index;
    public static bool IsValuedField => true;

    public SPoly(int ord, char x, K k)
    {
        Ord = ord;
        Poly = new KPoly<K>(x, k);
        Hash = (Poly.Degree, order: ord).GetHashCode();
    }

    public SPoly(int ord, KPoly<K> poly)
    {
        Ord = ord;
        Poly = new KPoly<K>(x, poly.KZero, poly.Coefs.Take(Ord).TrimSeq().ToArray());
        Hash = (Poly.Degree, order: ord).GetHashCode();
    }

    public int GetHashCodeSlow() => (Order: Ord, Poly.GetHashCodeSlow()).GetHashCode();
    public bool Equals(SPoly<K> other) => Poly.Equals(other.Poly); // Avoid collisions

    public int CompareTo(SPoly<K> other) => Poly.CompareTo(other.Poly);

    public int P => Poly.P;
    public SPoly<K> Clone => new(Ord, Poly);

    public SPoly<K> Inv()
    {
        if (Poly[0].IsZero())
            throw new DivideByZeroException();
        
        return new(Ord, Ring.NewtonInverse(Poly, Ord));
    }

    public IEnumerator<K> GetEnumerator() => Poly.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    public bool Invertible() => !Poly[0].IsZero();

    public SPoly<K> KMul(K k) => new(Ord, Poly.KMul(k));

    public int Hash { get; }
    public char x => Poly.x;
    public bool IsZero() => Poly.IsZero() || Poly[Ord].IsOne();
    public K this[int idx] => Poly[idx];
    public int Degree => Poly.Degree;

    public K KZero => Poly.KZero;
    public K KOne => Poly.KOne;
    public SPoly<K> Zero => new(Ord, Poly.Zero);
    public SPoly<K> One => new(Ord, Poly.One);
    public SPoly<K> X => new(Ord, Poly.X);
    public SPoly<K> Derivative => new(Ord, Poly.Derivative);
    public SPoly<K> Substitute(KPoly<K> f) => new(Ord, Poly.Substitute(f));
    public SPoly<K> Substitute(SPoly<K> f) => new(Ord, Poly.Substitute(f.Poly));

    public FracPoly<K> Substitute(FracPoly<K> f)
    {
        return Poly.Coefs.Select((k, i) => k * f.Pow(i)).Aggregate(f.Zero, (acc, a) => acc + a);
    }

    public SPoly<K> Compose(int k)
    {
        if (k < 0)
            throw new();

        if (k == 0)
            return X;

        return Substitute(Compose(k - 1));
    }

    public SPoly<K> Add(SPoly<K> e)
    {
        return new(Ord, Poly.Add(e.Poly));
    }

    public SPoly<K> Sub(SPoly<K> e)
    {
        return new(Ord, Poly.Sub(e.Poly));
    }

    public SPoly<K> Opp() => new(Ord, Poly.Opp());

    public SPoly<K> Mul(SPoly<K> e) => new(Ord, Poly.Mul(e.Poly));

    public (SPoly<K> quo, SPoly<K> rem) Div(SPoly<K> e)
    {
        var gcd = Ring.FastGCD(Poly, e.Poly);
        var e0 = new SPoly<K>(Ord, Poly / gcd);
        var e1 = new SPoly<K>(Ord, e.Poly / gcd);
        return (e0.Mul(e1.Inv()), Zero);
    }

    public SPoly<K> Mul(int k) => new(Ord, Poly.Mul(k));

    public SPoly<K> Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        if (Degree == 0)
            return new(Ord, Poly.Pow(k));

        return Ring.FastPow(this, k);
    }

    public KPoly<SPoly<K>> ToKPoly(char c) => new(c, this);
    public FracPoly<SPoly<K>> ToFracPoly(char c) => new(ToKPoly(c));

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsZero())
            return "0";

        return $"{Poly}";
    }

    public static SPoly<K> operator +(SPoly<K> a, SPoly<K> b) => a.Add(b);
    public static SPoly<K> operator +(int a, SPoly<K> b) => b.Add(b.One.Mul(a));
    public static SPoly<K> operator +(SPoly<K> a, int b) => a.Add(a.One.Mul(b));
    public static SPoly<K> operator -(SPoly<K> a) => a.Opp();
    public static SPoly<K> operator -(SPoly<K> a, SPoly<K> b) => a + (-b);
    public static SPoly<K> operator -(int a, SPoly<K> b) => a + (-b);
    public static SPoly<K> operator -(SPoly<K> a, int b) => a + (-b);
    public static SPoly<K> operator *(SPoly<K> a, SPoly<K> b) => a.Mul(b);
    public static SPoly<K> operator *(SPoly<K> a, int b) => a.Mul(b);
    public static SPoly<K> operator *(int a, SPoly<K> b) => b.Mul(a);
    public static SPoly<K> operator /(SPoly<K> a, SPoly<K> b) => a.Div(b).quo;
    public static SPoly<K> operator /(SPoly<K> a, int b) => a.Div(a.One.Mul(b)).quo;
    public static SPoly<K> operator /(int a, SPoly<K> b) => b.Inv().Mul(a);

    public static SPoly<K> operator +(SPoly<K> a, K b) => a + a.One.KMul(b);
    public static SPoly<K> operator +(K a, SPoly<K> b) => b.One.KMul(a) + b;
    public static SPoly<K> operator -(SPoly<K> a, K b) => a - a.One.KMul(b);
    public static SPoly<K> operator -(K a, SPoly<K> b) => b.One.KMul(a) - b;
    public static SPoly<K> operator *(SPoly<K> a, K b) => a.KMul(b);
    public static SPoly<K> operator *(K a, SPoly<K> b) => b.KMul(a);
    public static SPoly<K> operator /(SPoly<K> a, K b) => a.KMul(b.Inv());
    public static SPoly<K> operator /(K a, SPoly<K> b) => b.Inv().KMul(a);
    public static SPoly<K> operator +(SPoly<K> a, KPoly<K> b) => a + new SPoly<K>(a.Ord, b);
    public static SPoly<K> operator +(KPoly<K> a, SPoly<K> b) => new SPoly<K>(b.Ord, a) + b;
    public static SPoly<K> operator -(SPoly<K> a, KPoly<K> b) => a - new SPoly<K>(a.Ord, b);
    public static SPoly<K> operator -(KPoly<K> a, SPoly<K> b) => new SPoly<K>(b.Ord, a) - b;
    public static SPoly<K> operator *(SPoly<K> a, KPoly<K> b) => a * new SPoly<K>(a.Ord, b);
    public static SPoly<K> operator *(KPoly<K> a, SPoly<K> b) => new SPoly<K>(b.Ord, a) * b;
    public static SPoly<K> operator /(SPoly<K> a, KPoly<K> b) => a.Mul(new SPoly<K>(a.Ord, b).Inv());
    public static SPoly<K> operator /(KPoly<K> a, SPoly<K> b) => new SPoly<K>(b.Ord, a) * b.Inv();
}