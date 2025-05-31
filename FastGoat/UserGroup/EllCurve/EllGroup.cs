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

        Field = typeof(T).Name;
        if (Disc is Rational)
            Field = "Q";
        else if (Disc is ZnInt _a1)
            Field = $"Z/{_a1.P}Z";
        else if (Disc is ZnBigInt _a2)
            Field = $"Z/{_a2.Mod}Z";
        else if (Disc is ZnInt disc1)
            Field = $"Z/{disc1.P}Z";
        else if (Disc is ZnBigInt disc2)
            Field = $"Z/{disc2.Mod}Z";
        else if (Disc is EPoly<ZnInt> disc3)
        {
            if (disc3.F.Degree > 1)
                Field = $"GF({disc3.P}^{disc3.F.Degree})";
            else
                Field = $"GF({disc3.P})";
        }

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
            var (lhs, rhs) = Eq();
            return $"Elliptic curve {lhs} = {rhs}";
        }
    }

    public (Polynomial<T, Xi> lhs, Polynomial<T, Xi> rhs) Eq(Indeterminates<Xi> ind, Xi x, Xi y)
    {
        var X = new Polynomial<T, Xi>(new Monom<Xi>(ind, x), a1.One);
        var Y = new Polynomial<T, Xi>(new Monom<Xi>(ind, y), a1.One);
        var lhs = Y * Y + a1 * X * Y + a3 * Y;
        var rhs = X.Pow(3) + a2 * X * X + a4 * X + a6;
        return (lhs, rhs);
    }

    public (Polynomial<T, Xi> lhs, Polynomial<T, Xi> rhs) Eq(Polynomial<T, Xi> x, Polynomial<T, Xi> y) =>
        Eq(x.Indeterminates, x.ExtractIndeterminate, y.ExtractIndeterminate);

    public (Polynomial<T, Xi> lhs, Polynomial<T, Xi> rhs) Eq()
    {
        var (y, x) = Ring.Polynomial(Disc, MonomOrder.Lex, "Y", "X").Deconstruct();
        return Eq(x, y);
    }
    public (Polynomial<T, Xi> Eq, Polynomial<T, Xi> sd) GetPolynomials(Polynomial<T, Xi> x, Polynomial<T, Xi> y)
    {
        var (lhs, rhs) = Eq(x, y);
        return (lhs - rhs, 2 * y + a1 * x + a3);
    }

    public (Polynomial<T, Xi> X, Polynomial<T, Xi> Y, Polynomial<T, Xi> Eq, Polynomial<T, Xi> sd) GetPolynomials()
    {
        var (y, x) = Ring.Polynomial(Disc, MonomOrder.Lex, "Y", "X").Deconstruct();
        var (lhs, rhs) = Eq(x, y);
        return (x, y, lhs - rhs, 2 * y + a1 * x + a3);
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

    public bool Contains(T X, T Y)
    {
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
        if (!x1.Equals(x2))
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) + a1 * alpha - a2 - x2 - x1;
            var y3 = x3 * alpha - x1 * alpha + y1;
            return new(x3, -a1 * x3 - a3 - y3);
        }
        else
        {
            if (!y1.Equals(-a1 * x2 - a3 - y2))
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

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}