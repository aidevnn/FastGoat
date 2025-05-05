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
        var (A1, A2, A3, A4, A6) = Coefs = ellCoefs.Model;
        Coefs = ellCoefs.Model;
        Disc = ellCoefs.Disc;

        if (Disc.P != 0 && Disc.P <= 3)
            throw new GroupException(GroupExceptionType.GroupDef);

        var longCoefs = ellCoefs.Flat().ToLongWeierstrassForm();
        var (_, A2l, _, A4l, A6l) = longCoefs.Model;
        var (rl, sl, tl, ul) = longCoefs.TransCoef;
        LongForm = (A2l, A4l, A6l, rl, sl, tl, ul);

        var shortCoefs = ellCoefs.Flat().ToShortWeierstrassForm();
        var (_, _, _, A4s, A6s) = shortCoefs.Model;
        var (rs, ss, ts, us) = shortCoefs.TransCoef;
        ShortForm = (A4s, A6s, rs, ss, ts, us);

        Field = typeof(T).Name;
        if (A1 is Rational)
        if (Disc is Rational)
            Field = "Q";
        else if (A1 is ZnInt _a1)
            Field = $"Z/{_a1.P}Z";
        else if (A1 is ZnBigInt _a2)
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

        Hash = (Coefs, ShortForm).GetHashCode();
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

    public EllPt<T> ConvertToShort(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return new(x, y);
    }

    public EllPt<T> ConvertFromShort(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return new(x, y);
    }

    public EllPt<T> ConvertToLong(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return new(x, y);
    }

    public EllPt<T> ConvertFromLong(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return new(x, y);
    }

    public int Hash { get; }

    public string Eq
    {
        get
        {
            var (x, y) = Ring.Polynomial(Disc, "x", "y").Deconstruct();
            var lhs = y * y + a1 * x * y + a3 * y;
            var rhs = x.Pow(3) + a2 * x * x + a4 * x + a5;
            return $"Elliptic curve {lhs} = {rhs}";
        }
    }

    public string EqShort
    {
        get
        {
            var (x, y) = Ring.Polynomial(Disc, "x", "y").Deconstruct();
            var lhs = y * y;
            var rhs = x.Pow(3) + ShortForm.A * x + ShortForm.B;
            return $"Elliptic curve short form {lhs} = {rhs}";
        }
    }

    public string EqLong
    {
        get
        {
            var (x, y) = Ring.Polynomial(Disc, "x", "y").Deconstruct();
            var lhs = y * y;
            var rhs = x.Pow(3) + LongForm.A * x * x + LongForm.B * x + LongForm.C;
            return $"Elliptic curve long form {lhs} = {rhs}";
        }
    }

    public string Name => a1.IsZero() && a2.IsZero() && a3.IsZero()
        ? $"Ell[{a4},{a5}]({Field})".Replace(" ", "")
        : $"Ell[{a1},{a2},{a3},{a4},{a5}]({Field})".Replace(" ", "");

    public string Field { get; set; }
    public T Disc { get; }
    public (T a1, T a2, T a3, T a4, T a5) Coefs { get; }
    public (T A, T B, T r, T s, T t, T u) ShortForm { get; }
    public (T A, T B, T C, T r, T s, T t, T u) LongForm { get; }
    public string ShortFormStr => $"Ell[{ShortForm.A},{ShortForm.B}]({Field})".Replace(" ", "");
    public string LongFormStr => $"Ell[0,{LongForm.A},0,{LongForm.B},{LongForm.C}]({Field})".Replace(" ", "");
    public T a1 => Coefs.a1;
    public T a2 => Coefs.a2;
    public T a3 => Coefs.a3;
    public T a4 => Coefs.a4;
    public T a5 => Coefs.a5;

    public EllPt<T> this[params ValueType[] us]
    {
        get
        {
            if (us.Length == 2 && us[0] is int x0 && us[1] is int y0)
            {
                var (x1, y1) = (ShortForm.A.One * x0, ShortForm.A.One * y0);
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
        var rhs = X.Pow(3) + a2 * X * X + a4 * X + a5;
        return lhs.Equals(rhs);
    }

    public EllPt<T> O => new();
    public EllPt<T> Neutral() => O;

    public EllPt<T> Invert(EllPt<T> P)
    {
        if (P.IsO)
            return P;

        if (!Contains(P.X, P.Y))
            throw new GroupException(GroupExceptionType.GroupDef);

        var P1 = ConvertToShort(P);
        return ConvertFromShort(new EllPt<T>(P1.X, -P1.Y));
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

        var (e1_, e2_) = (ConvertToShort(e1), ConvertToShort(e2));
        var (x1, y1, x2, y2) = (e1_.X, e1_.Y, e2_.X, e2_.Y);
        if (!x1.Equals(x2))
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) - x1 - x2;
            var y3 = -y1 + alpha * (x1 - x3);
            return ConvertFromShort(new(x3, y3));
        }
        else
        {
            if (!y1.Equals(y2.Opp()))
            {
                var alpha = (3 * x1.Pow(2) + ShortForm.A) / (2 * y1);
                var x3 = alpha.Pow(2) - 2 * x1;
                var y3 = -y1 + alpha * (x1 - x3);
                return ConvertFromShort(new(x3, y3));
            }
            else
                return new();
        }
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}