using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup.EllCurve;

public readonly struct EllPoly<K> : IElt<EllPoly<K>>, IRingElt<EllPoly<K>>, IFieldElt<EllPoly<K>>,
    IModuleElt<K, EllPoly<K>>, IVsElt<K, EllPoly<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public Polynomial<K, Xi> Num { get; }
    public Polynomial<K, Xi> Denom { get; }
    public (Polynomial<K, Xi> eqEll, Polynomial<K, Xi> sd, Polynomial<K, Xi> dvp) Reduction { get; }
    public Polynomial<K, Xi> EqEll => Reduction.eqEll;
    public Polynomial<K, Xi> DivPol => Reduction.dvp;
    public Polynomial<K, Xi> SD => Reduction.sd;
    public bool IsNumber => !IsInfty && !IsIndeterminate;

    public EllPoly((Polynomial<K, Xi> eqEll, Polynomial<K, Xi> sd, Polynomial<K, Xi> dvp) red, Polynomial<K, Xi> num,
        Polynomial<K, Xi> denom)
    {
        (Reduction, Num, Denom) = (red, num, denom);

        if (denom.IsZero() && !num.IsZero())
            (Num, Denom) = (num.One, denom);
        else if (denom.IsZero())
            (Num, Denom) = (num.One, denom);
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

    public bool Equals(EllPoly<K> other)
    {
        return Sub(other).IsZero();
    }

    public int CompareTo(EllPoly<K> other)
    {
        if (IsIndeterminate)
            return -1;
        if (other.IsIndeterminate)
            return 1;
        if (IsInfty)
            return 1;
        if (other.IsInfty)
            return -1;
        return (Num * other.Denom).CompareTo(Denom * other.Num);
    }

    public bool IsInfty => !Num.IsZero() && Denom.IsZero();
    public bool IsIndeterminate => Num.IsZero() && Denom.IsZero();
    public EllPoly<K> Infty => new(Reduction, Num.One, Num.Zero);
    public EllPoly<K> Indeterminate => new(Reduction, Num.Zero, Num.Zero);
    public int Hash { get; }
    public bool IsZero() => !Denom.IsZero() && Num.IsZero();

    public bool IsDivZero()
    {
        if (Num.IsZero())
            return true;

        var num = Num.Div(DivPol).rem;
        var (y, x) = num.Indeterminates.Deconstruct();
        var dvp = DivPol.ToKPoly(x);
        var num0 = Ring.Decompose(num, y).Item1.ToDictionary(e => e.Key,
                e => Ring.FastGCD(e.Value.ToKPoly(x), dvp).Degree > 0 ? e.Value.Zero : e.Value)
            .Select(e => e.Key * e.Value).ToVec().Sum();

        return num0.IsZero();
    }

    public EllPoly<K> Zero => new(Reduction, Num.Zero, Num.One);
    public EllPoly<K> One => new(Reduction, Num.One, Num.One);

    public K KZero => Num.KZero;
    public K KOne => Num.KOne;

    public EllPoly<K> KMul(K k)
    {
        if (IsIndeterminate)
            return this;
        if (IsInfty)
            return k.IsZero() ? Indeterminate : Infty;

        return new(Reduction, k * Num, Denom);
    }

    public EllPoly<K> Add(EllPoly<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return Indeterminate;

        if (IsNumber && e.IsNumber)
        {
            var num = (Num * e.Denom + Denom * e.Num);
            var denom = Denom.Mul(e.Denom);
            return Simplify(Reduction, num, denom);
        }

        return Infty;
    }

    public EllPoly<K> Sub(EllPoly<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return Indeterminate;

        if (IsNumber && e.IsNumber)
        {
            var num = (Num * e.Denom - Denom * e.Num);
            var denom = Denom.Mul(e.Denom);
            return Simplify(Reduction, num, denom);
        }

        return Infty;
    }

    public EllPoly<K> Opp() => new(Reduction, -Num, Denom);

    public EllPoly<K> Mul(EllPoly<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return Indeterminate;

        if (IsNumber && e.IsNumber)
        {
            var num = (Num * e.Num);
            var denom = (Denom * e.Denom);
            return Simplify(Reduction, num, denom);
        }

        if (IsInfty && e.IsZero())
            return Indeterminate;

        if (e.IsInfty && IsZero())
            return Indeterminate;

        return Infty;
    }

    public (EllPoly<K> quo, EllPoly<K> rem) Div(EllPoly<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return (Indeterminate, Indeterminate);

        if (IsNumber && e.IsNumber)
        {
            var num = Num * e.Denom;
            var denom = Denom * e.Num;
            return (Simplify(Reduction, num, denom), Zero);
        }

        if (IsInfty && e.IsInfty)
            return (Indeterminate, Indeterminate);

        if (IsNumber && e.IsInfty)
            return (Zero, Zero);

        return (Infty, Infty);
    }

    public EllPoly<K> Mul(int k) => KMul(k * KOne);

    public EllPoly<K> Pow(int k)
    {
        if (IsIndeterminate)
            return Indeterminate;

        if (IsInfty)
            return k == 0 ? Indeterminate : k < 0 ? Zero : Infty;

        if (k == 0)
            return One;

        if (k < 0)
            return Inv().Pow(-k);

        return Ring.FastPow(this, k);
    }

    public int P => Num.P;

    public EllPoly<K> Inv()
    {
        if (IsNumber)
            return new(Reduction, Denom, Num);

        if (IsIndeterminate)
            return Indeterminate;

        return new(Reduction, Num.Zero, Num.One);
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsIndeterminate)
            return "?ind?";
        if (IsInfty)
            return "infty";

        if (Denom.Equals(Denom.One))
            return Num.ToString();

        return $"({Num})/({Denom})";
    }

    public FracPoly<FracPoly<K>> ToFracPoly(Xi x, Xi y)
    {
        var num = Num.ToFracPoly(x, y);
        var denom = Denom.ToFracPoly(x, y);
        return new(num, denom);
    }

    public bool Invertible() => !IsIndeterminate;

    public static (EllPoly<K> X, EllPoly<K> Y) GetXY(Polynomial<K, Xi> eqEll, Polynomial<K, Xi> sd,
        Polynomial<K, Xi> dvp)
    {
        var (Y, X) = eqEll.AllVariables
            .Select(xi => new EllPoly<K>((eqEll, sd, dvp), xi, xi.One))
            .Deconstruct();
        return (X, Y);
    }

    public static EllPoly<K> Simplify((Polynomial<K, Xi> eqEll, Polynomial<K, Xi> sd, Polynomial<K, Xi> dvp) red,
        Polynomial<K, Xi> num, Polynomial<K, Xi> denom)
    {
        if (denom.IsZero())
            return new(red, num, denom);

        var (y, x) = num.Indeterminates.Deconstruct();

        if (denom.DegreeOf(y) % 2 == 1)
        {
            num *= red.sd;
            denom *= red.sd;
        }

        num = num.Div(red.eqEll).rem;
        denom = denom.Div(red.eqEll).rem;
        if (!red.dvp.IsZero())
        {
            num = num.Div(red.dvp).rem;
            denom = denom.Div(red.dvp).rem;
        }
        else
        {
            var decNum = Ring.Decompose(num, y).Item1.ToDictionary(e => e.Key, e => (e.Value, e.Value.ToKPoly(x)));
            var decDenom = Ring.Decompose(denom, y).Item1.ToDictionary(e => e.Key, e => (e.Value, e.Value.ToKPoly(x)));
            var arr = decNum.Values.Select(e => e.Item2).Concat(decDenom.Values.Select(e => e.Item2))
                .Where(e => !e.IsZero()).Distinct().ToArray();
            if (arr.Length != 0)
            {
                var gcd = Ring.FastGCD(arr).ToPolynomial(num.Indeterminates, x);
                num = decNum.Select(e => e.Key * (e.Value.Value / gcd)).ToVec().Sum();
                denom = decDenom.Select(e => e.Key * (e.Value.Value / gcd)).ToVec().Sum();
            }
        }

        return new(red, num, denom);
    }

    public static EllPoly<K> operator +(EllPoly<K> a, EllPoly<K> b) => a.Add(b);

    public static EllPoly<K> operator +(int a, EllPoly<K> b) => (b.One.Mul(a)).Add(b);

    public static EllPoly<K> operator +(EllPoly<K> a, int b) => a.Add(a.One.Mul(b));

    public static EllPoly<K> operator -(EllPoly<K> a) => a.Opp();

    public static EllPoly<K> operator -(EllPoly<K> a, EllPoly<K> b) => a.Sub(b);

    public static EllPoly<K> operator -(int a, EllPoly<K> b) => (b.One.Mul(a)).Sub(b);

    public static EllPoly<K> operator -(EllPoly<K> a, int b) => a.Sub(a.One.Mul(b));

    public static EllPoly<K> operator *(EllPoly<K> a, EllPoly<K> b) => a.Mul(b);

    public static EllPoly<K> operator *(int a, EllPoly<K> b) => b.Mul(a);

    public static EllPoly<K> operator *(EllPoly<K> a, int b) => a.Mul(b);

    public static EllPoly<K> operator /(EllPoly<K> a, EllPoly<K> b) => a.Div(b).quo;

    public static EllPoly<K> operator /(EllPoly<K> a, int b) => new(a.Reduction, a.Num, a.Denom * b);

    public static EllPoly<K> operator /(int a, EllPoly<K> b) => a * b.Inv();

    public static EllPoly<K> operator +(EllPoly<K> a, K b) => a + a.One.KMul(b);

    public static EllPoly<K> operator +(K a, EllPoly<K> b) => b.One.KMul(a) + b;

    public static EllPoly<K> operator -(EllPoly<K> a, K b) => a - a.One.KMul(b);

    public static EllPoly<K> operator -(K a, EllPoly<K> b) => b.One.KMul(a) - b;

    public static EllPoly<K> operator *(EllPoly<K> a, K b) => a.KMul(b);

    public static EllPoly<K> operator *(K a, EllPoly<K> b) => b.KMul(a);

    public static EllPoly<K> operator /(EllPoly<K> a, K b) => new(a.Reduction, a.Num, a.Denom * b);
    public static double Abs(EllPoly<K> t) => Polynomial<K, Xi>.Abs(t.Num) / Polynomial<K, Xi>.Abs(t.Denom);

    public static bool IsValuedField => K.IsValuedField;
}

public static class EllPolyExt
{
    public static EllPt<FracPoly<FracPoly<K>>> ToFracPoly<K>(this EllPt<EllPoly<K>> P, Xi x, Xi y)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (P.IsO)
            return new();

        return new(P.X.ToFracPoly(x, y), P.Y.ToFracPoly(x, y));
    }

    public static EllPoly<K> ToEllPoly<K>(this FracPoly<FracPoly<K>> f, EllPoly<K> x, EllPoly<K> y)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var num = f.Num.Coefs.Select((c, i) => c.Num.Substitute(x.Num) * y.Num.Pow(i)).ToVec().Sum();
        var denom = f.Num.Coefs.Select(c => c.Denom.Substitute(x.Num)).ToVec().Sum();
        return EllPoly<K>.Simplify(x.Reduction, num, denom);
    }

    public static EllPt<EllPoly<K>> ToEllPoly<K>(this EllPt<FracPoly<FracPoly<K>>> P, EllPoly<K> x, EllPoly<K> y)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (P.IsO)
            return new();

        return new(P.X.ToEllPoly(x, y), P.Y.ToEllPoly(x, y));
    }
}