using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.EllCurve;

public readonly struct EllFracPoly<K> : IElt<EllFracPoly<K>>, IRingElt<EllFracPoly<K>>, IFieldElt<EllFracPoly<K>>,
    IModuleElt<K, EllFracPoly<K>>, IVsElt<K, EllFracPoly<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public EllPoly<K> Num { get; }
    public EllPoly<K> Denom { get; }
    public (EllPoly<K> eqEll, EllPoly<K> sd, EllPoly<K> dvp) Reduction { get; }
    public EllPoly<K> EqEll => Reduction.eqEll;
    public EllPoly<K> DivPol => Reduction.dvp;
    public EllPoly<K> SD => Reduction.sd;

    public EllFracPoly((EllPoly<K> eqEll, EllPoly<K> sd, EllPoly<K> dvp) red, EllPoly<K> num,
        EllPoly<K> denom)
    {
        (Reduction, Num, Denom) = (red, num, denom);

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

        Hash = (Num.Hash, Denom.Hash, Reduction.GetHashCode()).GetHashCode();
    }

    public bool Equals(EllFracPoly<K> other)
    {
        return Sub(other).IsZero();
    }

    public int CompareTo(EllFracPoly<K> other)
    {
        return (Num * other.Denom).CompareTo(Denom * other.Num);
    }

    public int Hash { get; }
    public bool IsZero() => !Denom.IsZero() && Num.IsZero();

    public bool IsDivZero()
    {
        if (IsZero())
            return true;

        if (DivPol.Degree == 0)
            return false;

        var fdiv = DivPol;
        var arr = Num.DecomposeX2().Values.Where(e => !e.IsZero());
        if (arr.All(f => f.Degree > 0 && EC.FastGCD(f, fdiv).Degree > 0))
            return true;

        return false;
    }

    public EllFracPoly<K> Zero => new(Reduction, Num.Zero, Num.One);
    public EllFracPoly<K> One => new(Reduction, Num.One, Num.One);

    public K KZero => Num.KZero;
    public K KOne => Num.KOne;

    public EllFracPoly<K> KMul(K k)
    {
        return new(Reduction, k * Num, Denom);
    }

    public EllFracPoly<K> Add(EllFracPoly<K> e)
    {
        var num = (Num * e.Denom + Denom * e.Num);
        var denom = Denom.Mul(e.Denom);
        return Simplify(Reduction, num, denom);
    }

    public EllFracPoly<K> Sub(EllFracPoly<K> e)
    {
        var num = (Num * e.Denom - Denom * e.Num);
        var denom = Denom.Mul(e.Denom);
        return Simplify(Reduction, num, denom);
    }

    public EllFracPoly<K> Opp() => new(Reduction, -Num, Denom);

    public EllFracPoly<K> Mul(EllFracPoly<K> e)
    {
        var num = (Num * e.Num);
        var denom = (Denom * e.Denom);
        return Simplify(Reduction, num, denom);
    }

    public (EllFracPoly<K> quo, EllFracPoly<K> rem) Div(EllFracPoly<K> e)
    {
        var num = Num * e.Denom;
        var denom = Denom * e.Num;
        return (Simplify(Reduction, num, denom), Zero);
    }

    public EllFracPoly<K> Mul(int k) => KMul(k * KOne);

    public EllFracPoly<K> Pow(int k)
    {
        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        return this.FastPow(k);
    }

    public int P => Num.P;

    public EllFracPoly<K> Inv()
    {
        return new(Reduction, Denom, Num);
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (Denom.Equals(Denom.One))
            return Num.ToString();

        return $"({Num})/({Denom})";
    }

    public bool Invertible() => !IsDivZero();

    public static EllFracPoly<K> Simplify((EllPoly<K> eqEll, EllPoly<K> sd, EllPoly<K> dvp) red, EllPoly<K> num,
        EllPoly<K> denom)
    {
        if (denom.DegreeOfX2 % 2 == 1)
        {
            num *= red.sd;
            denom *= red.sd;
        }

        num = num.Div(red.eqEll).rem;
        denom = denom.Div(red.eqEll).rem;
        if (red.dvp.Degree > 0)
        {
            num = num.Div(red.dvp).rem;
            denom = denom.Div(red.dvp).rem;
        }
        else
        {
            var decNum = num.DecomposeX2();
            var decDenom = denom.DecomposeX2();
            var arr = decNum.Values.Concat(decDenom.Values).Where(e => e.Degree > 0).Distinct().ToArray();
            if (arr.Length > 1)
            {
                var gcd = EC.FastGCD(arr).Monic;
                num = decNum.Select(e => e.Key * (e.Value / gcd)).Aggregate(gcd.Zero, (acc, e) => e + acc);
                denom = decDenom.Select(e => e.Key * (e.Value / gcd)).Aggregate(gcd.Zero, (acc, e) => e + acc);
            }
        }

        return new(red, num, denom);
    }

    public static EllFracPoly<K> operator +(EllFracPoly<K> a, EllFracPoly<K> b) => a.Add(b);

    public static EllFracPoly<K> operator +(int a, EllFracPoly<K> b) => (b.One.Mul(a)).Add(b);

    public static EllFracPoly<K> operator +(EllFracPoly<K> a, int b) => a.Add(a.One.Mul(b));

    public static EllFracPoly<K> operator -(EllFracPoly<K> a) => a.Opp();

    public static EllFracPoly<K> operator -(EllFracPoly<K> a, EllFracPoly<K> b) => a.Sub(b);

    public static EllFracPoly<K> operator -(int a, EllFracPoly<K> b) => (b.One.Mul(a)).Sub(b);

    public static EllFracPoly<K> operator -(EllFracPoly<K> a, int b) => a.Sub(a.One.Mul(b));

    public static EllFracPoly<K> operator *(EllFracPoly<K> a, EllFracPoly<K> b) => a.Mul(b);

    public static EllFracPoly<K> operator *(int a, EllFracPoly<K> b) => b.Mul(a);

    public static EllFracPoly<K> operator *(EllFracPoly<K> a, int b) => a.Mul(b);

    public static EllFracPoly<K> operator /(EllFracPoly<K> a, EllFracPoly<K> b) => a.Div(b).quo;

    public static EllFracPoly<K> operator /(EllFracPoly<K> a, int b) => new(a.Reduction, a.Num, a.Denom * b);

    public static EllFracPoly<K> operator /(int a, EllFracPoly<K> b) => a * b.Inv();

    public static EllFracPoly<K> operator +(EllFracPoly<K> a, K b) => a + a.One.KMul(b);

    public static EllFracPoly<K> operator +(K a, EllFracPoly<K> b) => b.One.KMul(a) + b;

    public static EllFracPoly<K> operator -(EllFracPoly<K> a, K b) => a - a.One.KMul(b);

    public static EllFracPoly<K> operator -(K a, EllFracPoly<K> b) => b.One.KMul(a) - b;

    public static EllFracPoly<K> operator *(EllFracPoly<K> a, K b) => a.KMul(b);

    public static EllFracPoly<K> operator *(K a, EllFracPoly<K> b) => b.KMul(a);

    public static EllFracPoly<K> operator /(EllFracPoly<K> a, K b) => new(a.Reduction, a.Num, a.Denom * b);
    public static double Abs(EllFracPoly<K> t) => EllPoly<K>.Abs(t.Num) / EllPoly<K>.Abs(t.Denom);

    public static bool IsValuedField => false;
}