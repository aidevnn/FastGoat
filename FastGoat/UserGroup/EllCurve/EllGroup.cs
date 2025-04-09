using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.EllCurve;

// Daniel Guin - Thomas Hausberger, Algebre Tome 1, page 166

public struct EllGroup<T> : IGroup<EllPt<T>> where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>
{
    public EllGroup(T a1, T a2, T a3, T a4, T a5)
    {
        var A1 = -a1.Pow(4) / 48 - a1.Pow(2) * a2 / 6 + a1 * a3 / 2 - a2.Pow(2) / 3 + a4;
        var B1 = a1.Pow(6) / 864 + a1.Pow(4) * a2 / 72 - a1.Pow(3) * a3 / 24 + a1.Pow(2) * a2.Pow(2) / 18 -
            a1.Pow(2) * a4 / 12 - a1 * a3 * a2 / 6 + 2 * a2.Pow(3) / 27 + a3.Pow(2) / 4 - a2 * a4 / 3 + a5;

        Disc = -4 * A1.Pow(3) - 27 * B1.Pow(2);
        if (Disc.IsZero())
            throw new GroupException(GroupExceptionType.GroupDef);

        T d11 = a1.One, d12 = a1.One;
        if (A1 is Rational _A1 && B1 is Rational _B1)
        {
            var sqDivs864 = new[] { 1, 4, 9, 16, 36, 144 } // Square Divisidors of 864 
                .Select(div => (div, pow2: div * div, pow3: div * div * div)).ToArray();
            var (sqDiv, _, sqDivPow3) = sqDivs864.OrderBy(f => f.div)
                .First(div => div.pow2 % _A1.Denom == 0 && div.pow3 % _B1.Denom == 0);
            dynamic _d11 = new Rational(sqDiv);
            dynamic _d12 = new Rational(IntExt.SqrtBigInt(sqDivPow3));
            (d11, d12) = (_d11, _d12);
        }

        A1 *= d11.Pow(2);
        B1 *= d11.Pow(3);

        var (A2, B2, C2) = (a1.Pow(2) / 4 + a2, a1 * a3 / 2 + a4, a3.Pow(2) / 4 + a5);

        T d21 = a1.One, d22 = a1.One;
        if (A2 is Rational A0 && B2 is Rational B0 && C2 is Rational C0 &&
            (!A0.IsInteger() || !B0.IsInteger() || !C0.IsInteger()))
            (d21, d22) = (4 * a1.One, 8 * a1.One);

        A2 *= d21;
        B2 *= d21.Pow(2);
        C2 *= d21.Pow(3);

        Coefs = (a1, a2, a3, a4, a5);
        ShortForm = (A1, B1, d11, d12);
        LongForm = (A2, B2, C2, d21, d22);

        Field = typeof(T).Name;
        if (A1 is Rational)
            Field = "Q";
        else if (A1 is ZnInt _a1)
            Field = $"Z/{_a1.P}Z";
        else if (A1 is ZnBigInt _a2)
            Field = $"Z/{_a2.Mod}Z";

        Name = a1.IsZero() && a2.IsZero() && a3.IsZero()
            ? $"Ell[{a4},{a5}]({Field})".Replace(" ", "")
            : $"Ell[{a1},{a2},{a3},{a4},{a5}]({Field})".Replace(" ", "");
        Hash = (Coefs, ShortForm).GetHashCode();
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

        var _x = (pt.X + a1.Pow(2) / 12 + a2 / 3) * ShortForm.d1;
        var _y = (pt.Y + (a1 * pt.X + a3) / 2) * ShortForm.d2;
        return new(_x, _y);
    }

    public EllPt<T> ConvertFromShort(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        var _x = pt.X / ShortForm.d1 - a1.Pow(2) / 12 - a2 / 3;
        var _y = pt.Y / ShortForm.d2;
        return new(_x, _y - (a1 * _x + a3) / 2);
    }

    public EllPt<T> ConvertToLong(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        return new(pt.X * LongForm.d1, (pt.Y + (a1 * pt.X + a3) / 2) * LongForm.d2);
    }

    public EllPt<T> ConvertFromLong(EllPt<T> pt)
    {
        if (pt.IsO)
            return pt;

        return new(pt.X / LongForm.d1, (pt.Y / LongForm.d2 - (a1 * pt.X / LongForm.d1 + a3) / 2));
    }

    public int Hash { get; }
    public string Name { get; }
    public string Field { get; }
    public T Disc { get; }
    public (T a1, T a2, T a3, T a4, T a5) Coefs { get; }
    public (T A, T B, T d1, T d2) ShortForm { get; }
    public (T A, T B, T C, T d1, T d2) LongForm { get; }
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