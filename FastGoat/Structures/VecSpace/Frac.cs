using System.Reflection;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Structures.VecSpace;

public readonly struct Frac<K> : IElt<Frac<K>>, IRingElt<Frac<K>>, IFieldElt<Frac<K>>, IModuleElt<K, Frac<K>>,
    IVsElt<K, Frac<K>>, IFieldInfElt<Frac<K>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public KPoly<K> Num { get; }
    public KPoly<K> Denom { get; }
    public bool IsInfinity { get; }
    public bool IsIndeterminate { get; }

    public static bool IsValuedField => KPoly<K>.IsValuedField;
    public static double Abs(Frac<K> e) => KPoly<K>.Abs(e.Num) / KPoly<K>.Abs(e.Denom);

    public Frac(char x, K scalar)
    {
        IsIndeterminate = IsInfinity = false;
        var num = new KPoly<K>(x, scalar);
        Num = num.X;
        Denom = Num.One;
        Hash = (Num.Hash, Denom.Hash, IsIndeterminate, IsInfinity).GetHashCode();
    }

    public Frac(KPoly<K> num)
    {
        (Num, IsIndeterminate, IsInfinity) = num.Rec();
        Denom = Num.One;
        Hash = (Num.Hash, Denom.Hash, IsIndeterminate, IsInfinity).GetHashCode();
    }

    public Frac(KPoly<K> num, KPoly<K> denom)
    {
        IsIndeterminate = IsInfinity = false;
        if (num.P != denom.P || num.x != denom.x)
        {
            Console.WriteLine($"Char:{num.P}/{denom.P} symb:{num.x}/{denom.x}");
            throw new GroupException(GroupExceptionType.GroupDef);
        }
        
        var (_, isIndNum, isInfNum) = num.Rec();
        var (_, isIndDenom, isInfDenom) = denom.Rec();

        if (isIndNum || isIndDenom)
        {
            IsIndeterminate = true;
            IsInfinity = false;
            (Num, Denom) = (num.Zero, denom.Zero);
            return;
        }
        else if (isInfNum && !isInfDenom)
        {
            IsIndeterminate = false;
            IsInfinity = true;
            (Num, Denom) = (num.One, denom.Zero);
            return;
        }
        else if (isInfDenom && !isInfNum)
        {
            IsIndeterminate = false;
            IsInfinity = false;
            (Num, Denom) = (num.Zero, denom.One);
            return;
        }
        else if (isInfNum && isInfDenom)
        {
            IsIndeterminate = true;
            IsInfinity = false;
            (Num, Denom) = (num.Zero, denom.Zero);
            return;
        }

        if (denom.IsZero())
        {
            IsIndeterminate = num.IsZero();
            IsInfinity = !IsIndeterminate;
            (Num, Denom) = (num.One, denom);
            return;
        }

        var (q, r) = num.Div(denom);
        if (r.IsZero())
        {
            Num = q;
            Denom = q.One;
        }
        else
        {
            // Console.WriteLine(new{num, denom});
            // Console.WriteLine($"num => {num.Rec()}");
            // Console.WriteLine($"denom => {denom.Rec()}");
            var gcd = Ring.FastGCD(num, denom);
            // Console.WriteLine(new { gcd });
            var num0 = num.Div(gcd).quo;
            var denom0 = denom.Div(gcd).quo;
            var c = denom0.Coefs.Last();
            Num = num0 / c;
            Denom = denom0 / c;
        }

        Hash = (Num.Hash, Denom.Hash, IsIndeterminate, IsInfinity).GetHashCode();
    }

    public char x => Num.x;

    public bool Equals(Frac<K> other)
    {
        if (IsIndeterminate || other.IsIndeterminate)
            return false;

        if (IsInfinity && !other.IsInfinity)
            return false;
        if (other.IsInfinity && !IsInfinity)
            return false;
        if (IsInfinity && other.IsInfinity)
            return true;

        return Num.Equals(other.Num) && Denom.Equals(other.Denom);
    }

    public int CompareTo(Frac<K> other)
    {
        if (IsIndeterminate)
            return -1;
        if (other.IsIndeterminate)
            return 1;
        if (IsInfinity)
            return 1;
        if (other.IsInfinity)
            return -1;
        return (Num * other.Denom).CompareTo(Denom * other.Num);
    }

    public K KZero => Num.KZero;
    public K KOne => Num.KOne;
    public Frac<K> Indeterminate => new(Num.Zero, Num.Zero);
    public Frac<K> Infinity => new(Num.One, Num.Zero);
    public bool IsDeterminate => !IsInfinity && !IsIndeterminate;
    public Frac<K> KMul(K k)
    {
        if (IsIndeterminate)
            return this;
        if (IsInfinity)
            return k.IsZero() ? Indeterminate : Infinity;
        
        return new(k * Num, Denom);
    }
    public int Hash { get; }

    public bool IsZero() => IsDeterminate && Num.IsZero();

    public int P => Num.P;
    public Frac<K> Zero => new(Num.Zero);
    public Frac<K> One => new(Num.One);
    public Frac<K> X => new(Num.X);

    public Frac<K> Add(Frac<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return Indeterminate;
        
        if (IsDeterminate && e.IsDeterminate)
            return new(Num.Mul(e.Denom).Add(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

        return Infinity;
    }

    public Frac<K> Sub(Frac<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return Indeterminate;
        
        if (IsDeterminate && e.IsDeterminate)
            return new(Num.Mul(e.Denom).Sub(e.Num.Mul(Denom)), Denom.Mul(e.Denom));

        return Infinity;
    }

    public Frac<K> Opp() => new(Num.Opp(), Denom);

    public Frac<K> Mul(Frac<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return Indeterminate;
        
        if (IsDeterminate && e.IsDeterminate)
            return new(Num * e.Num, Denom * e.Denom);

        if (IsInfinity && e.IsZero())
            return Indeterminate;

        if (e.IsInfinity && IsZero())
            return Indeterminate;

        return Infinity;
    }

    public Frac<K> Inv()
    {
        if (IsDeterminate)
            return new(Denom, Num);

        if (IsIndeterminate)
            return Indeterminate;

        return new(Num.Zero);
    }

    public bool Invertible() => IsDeterminate && !Num.IsZero();

    public (Frac<K> quo, Frac<K> rem) Div(Frac<K> e)
    {
        if (IsIndeterminate || e.IsIndeterminate)
            return (Indeterminate, Indeterminate);

        if (IsDeterminate && e.IsDeterminate)
            return (new Frac<K>(Num * e.Denom, Denom * e.Num), Zero);

        if (IsInfinity && e.IsInfinity)
            return (Indeterminate, Indeterminate);

        if (IsDeterminate && e.IsInfinity)
            return (Zero, Zero);

        return (Infinity, Infinity);
    }

    public Frac<K> Mul(int k) => KMul(k * KOne);

    public Frac<K> Pow(int k)
    {
        if (IsIndeterminate)
            return Indeterminate;

        if (IsInfinity)
            return k == 0 ? Indeterminate : k < 0 ? Zero : Infinity;
        
        if (k < 0)
            return new(Denom.Pow(-k), Num.Pow(-k));

        return new(Num.Pow(k), Denom.Pow(k));
    }

    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        if (IsIndeterminate)
            return "?ind?";
        if (IsInfinity)
            return "infty";
        
        if (Denom.Equals(Denom.One))
            return Num.ToString();

        return $"({Num})/({Denom})";
    }

    public static Frac<K> operator +(Frac<K> a, Frac<K> b) => a.Add(b);

    public static Frac<K> operator +(int a, Frac<K> b) => (b.One.Mul(a)).Add(b);

    public static Frac<K> operator +(Frac<K> a, int b) => a.Add(a.One.Mul(b));

    public static Frac<K> operator -(Frac<K> a) => a.Opp();

    public static Frac<K> operator -(Frac<K> a, Frac<K> b) => a.Sub(b);

    public static Frac<K> operator -(int a, Frac<K> b) => (b.One.Mul(a)).Sub(b);

    public static Frac<K> operator -(Frac<K> a, int b) => a.Sub(a.One.Mul(b));

    public static Frac<K> operator *(Frac<K> a, Frac<K> b) => a.Mul(b);

    public static Frac<K> operator *(int a, Frac<K> b) => b.Mul(a);

    public static Frac<K> operator *(Frac<K> a, int b) => a.Mul(b);

    public static Frac<K> operator /(Frac<K> a, Frac<K> b) => a.Div(b).quo;

    public static Frac<K> operator /(Frac<K> a, int b) => new(a.Num, a.Denom * b);

    public static Frac<K> operator /(int a, Frac<K> b) => new(b.Denom, a * b.Num);

    public static Frac<K> operator +(Frac<K> a, K b) => a + new Frac<K>(new KPoly<K>(a.x, b));

    public static Frac<K> operator +(K a, Frac<K> b) => new Frac<K>(new KPoly<K>(b.x, a)) + b;

    public static Frac<K> operator -(Frac<K> a, K b) => a - new Frac<K>(new KPoly<K>(a.x, b));

    public static Frac<K> operator -(K a, Frac<K> b) => new Frac<K>(new KPoly<K>(b.x, a)) - b;

    public static Frac<K> operator *(Frac<K> a, K b) => new(a.Num * b, a.Denom);

    public static Frac<K> operator *(K a, Frac<K> b) => new(a * b.Num, b.Denom);

    public static Frac<K> operator /(Frac<K> a, K b) => new(a.Num, a.Denom * b);

    public static Frac<K> operator +(Frac<K> a, KPoly<K> b) => a + new Frac<K>(b);

    public static Frac<K> operator +(KPoly<K> a, Frac<K> b) => new Frac<K>(a) + b;

    public static Frac<K> operator -(Frac<K> a, KPoly<K> b) => a - new Frac<K>(b);

    public static Frac<K> operator -(KPoly<K> a, Frac<K> b) => new Frac<K>(a) - b;

    public static Frac<K> operator *(Frac<K> a, KPoly<K> b) => new(a.Num * b, a.Denom);

    public static Frac<K> operator *(KPoly<K> a, Frac<K> b) => new(a * b.Num, b.Denom);

    public static Frac<K> operator /(Frac<K> a, KPoly<K> b) => new(a.Num, a.Denom * b);
}

file static class FracRec
{
    private static (KPoly<K>, bool isInd, bool isInf) Rec<K, T>(this KPoly<K> F) 
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
        where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
    {
        if (F.KOne is Frac<T> e)
        {
            if (F.Coefs.Cast<Frac<T>>().Any(f => f.IsIndeterminate))
                return ((dynamic)F.Zero * e.Indeterminate, true, false);
            
            if (F.Coefs.Cast<Frac<T>>().Any(f => f.IsInfinity))
                return ((dynamic)F.One * e.Infinity, false, true);
        }

        return (F, false, false);
    }

    public static (KPoly<K>, bool isInd, bool isInf) Rec<K>(this KPoly<K> F)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        if (F.KOne is Frac<Rational>)
            return Rec<K, Rational>(F);
        if (F.KOne is Frac<ZnInt>)
            return Rec<K, ZnInt>(F);
        if (F.KOne is Frac<EPoly<Rational>>)
            return Rec<K, EPoly<Rational>>(F);
        if (F.KOne is Frac<EPoly<ZnInt>>)
            return Rec<K, EPoly<ZnInt>>(F);

        return (F, false, false);
    }
}