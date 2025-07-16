using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.EllCurve;

// Daniel Guin - Thomas Hausberger, Algebre Tome 1, page 166

public struct EllGroup<T> : IGroup<EllPt<T>> where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
{
    public EllGroup(EllCoefs<T> ellCoefs)
    {
        Coefs = ellCoefs.Model;
        Coefs = ellCoefs.Model;
        Disc = ellCoefs.Disc;

        if (Disc.IsZero())
            throw new GroupException(GroupExceptionType.GroupDef);

        Field = GetFieldName(Disc);
        Hash = (Coefs, Field).GetHashCode();
    }

    public EllGroup(T a1, T a2, T a3, T a4, T a6) : this(new EllCoefs<T>(a1, a2, a3, a4, a6))
    {
    }

    public EllGroup(T a, T b) : this(a.Zero, a.Zero, a.Zero, a, b)
    {
    }

    public EllGroup(T a, T b, T c) : this(a.Zero, a, a.Zero, b, c)
    {
    }

    public IEnumerator<EllPt<T>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EllPt<T>>? other) => other?.Hash == Hash;

    public int Hash { get; }

    public string EqStr
    {
        get
        {
            var (_, Y, X) = EC.EllPoly(a1);
            var lhs = Y * Y + a1 * X * Y + a3 * Y;
            var rhs = X.Pow(3) + a2 * X * X + a4 * X + a6;
            return $"Elliptic curve {lhs} = {rhs}";
        }
    }

    public string Name => a1.IsZero() && a2.IsZero() && a3.IsZero()
        ? $"Ell[{a4},{a6}]({Field})".Replace(" ", "")
        : $"Ell[{a1},{a2},{a3},{a4},{a6}]({Field})".Replace(" ", "");

    public string Field { get; set; }
    public T Disc { get; }
    public (T a1, T a2, T a3, T a4, T a6) Coefs { get; }
    public T[] ArrCoefs => [a1, a2, a3, a4, a6];
    public T a1 => Coefs.a1;
    public T a2 => Coefs.a2;
    public T a3 => Coefs.a3;
    public T a4 => Coefs.a4;
    public T a6 => Coefs.a6;
    
    public bool CheckValidity { get; set; } = true;

    public EllPt<T> this[params ValueType[] us]
    {
        get
        {
            if (us.Length == 2 && us[0] is int x0 && us[1] is int y0)
            {
                var (x1, y1) = (a1.One * x0, a1.One * y0);
                if (!Contains(x1, y1))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<T>(x1, y1);
            }

            if (us.Length == 2 && us[0] is T x && us[1] is T y)
            {
                if (!Contains(x, y))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<T>(x, y);
            }

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EllPt<T>> GetElements()
    {
        yield return new EllPt<T>();
    }

    public IEnumerable<EllPt<T>> GetGenerators()
    {
        yield return new EllPt<T>();
    }

    public K Lhs<K>(K X, K Y) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>,
        IVsElt<T, K>, IModuleElt<T, K>
    {
        return Y * Y + a1 * X * Y + a3 * Y;
    }

    public K Rhs<K>(K X) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>,
        IVsElt<T, K>, IModuleElt<T, K>
    {
        return X.Pow(3) + a2 * X * X + a4 * X + a6;
    }

    public (K lhs, K rhs) Eq<K>(K X, K Y) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>,
        IVsElt<T, K>, IModuleElt<T, K>
    {
        return (Lhs(X, Y), Rhs(X));
    }
    
    public bool Contains(T X, T Y)
    {
        if (!CheckValidity)
            return true;

        var lhs = Y * Y + a1 * X * Y + a3 * Y;
        var rhs = X.Pow(3) + a2 * X * X + a4 * X + a6;
        return lhs.Equals(rhs);
    }

    public EllPt<T> O => new();
    public EllPt<T> Neutral() => O;

    public EllPt<T> Invert(EllPt<T> P)
    {
        if (P.IsO)
            return P;

        if (!Contains(P.X, P.Y))
        {
            Console.WriteLine(new { P, E = this });
            throw new GroupException(GroupExceptionType.GroupDef);
        }

        return new(P.X, -a1 * P.X - a3 - P.Y);
    }

    public EllPt<T> Op(EllPt<T> e1, EllPt<T> e2)
    {
        if (e1.IsO)
            return e2;

        if (e2.IsO)
            return e1;

        if (!Contains(e1.X, e1.Y) || !Contains(e2.X, e2.Y))
        {
            Console.WriteLine(new { e1, e2, E = this });
            throw new GroupException(GroupExceptionType.GroupDef);
        }

        var ((x1, y1), (x2, y2)) = (e1, e2);
        if ((x1 - x2).Invertible())
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) + a1 * alpha - a2 - x2 - x1;
            var y3 = x3 * alpha - x1 * alpha + y1;
            return new(x3, -a1 * x3 - a3 - y3);
        }
        else
        {
            if ((y1 + a1 * x2 + a3 + y2).Invertible() && (x1 * a1 + a3 + 2 * y1).Invertible())
            {
                var alpha = (3 * x1.Pow(2) - y1 * a1 + 2 * x1 * a2 + a4) / (x1 * a1 + a3 + 2 * y1);
                var x3 = alpha.Pow(2) + a1 * alpha - a2 - 2 * x1;
                var y3 = alpha * (x3 - x1) + y1;
                return new(x3, -a1 * x3 - a3 - y3);
            }
            else
                return new();
        }
    }

    public EllCoefs<T> ToEllCoefs() => new(a1, a2, a3, a4, a6);
    public EllCoefs<T> ToLongWeierstrassForm() => ToEllCoefs().ToLongWeierstrassForm();
    public EllCoefs<T> ToShortWeierstrassForm() => ToEllCoefs().ToShortWeierstrassForm();

    public EllGroup<K> ToEllGroup<K>(K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>,
        IVsElt<T, K>, IModuleElt<T, K>
    {
        var o = scalar.One;
        return new(a1 * o, a2 * o, a3 * o, a4 * o, a6 * o);
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;

    public static string GetFieldName(T e)
    {
        if (e is Rational)
            return "Q";
        else if (e is ZnInt e1)
            return $"Z/{e1.P}Z";
        else if (e is ZnBigInt e2)
            return $"Z/{e2.Mod}Z";
        else if (e is EPoly<ZnInt> e3)
        {
            if (e3.F.Degree > 1)
                return $"GF({e3.P}^{e3.F.Degree})";
            else
                return $"GF({e3.P})";
        }
        else if (e is TriVarFrac<Rational>)
            return "Q[X,Y]";
        else if (e is TriVarFrac<ZnInt> e4)
            return $"F{e4.KOne.Mod}[X,Y]";
        else if (e is TriVarFrac<ZnBigInt> e5)
            return $"F{e5.KOne.Mod}[X,Y]";
        else if (e is EPoly<Rational> e6)
            return $"Q({e6.F.x})";

        return typeof(T).Name;
    }
}