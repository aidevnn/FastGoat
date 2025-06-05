using System.Collections;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.EllCurve;


public struct EllGroupSymb<K> : IGroup<EllPt<EllFracPoly<K>>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public EllGroupSymb(EllCoefs<K> ellCoefs, EllPoly<K> divPol)
    {
        (Y, X) = ellCoefs.EllFracPolyYX(divPol);
        DivPolynomial = divPol;
        var ec = ellCoefs.ToEllCoefs(X);
        Coefs = ec.Model;
        Disc = ec.Disc;

        var discP = ellCoefs.Disc is ZnBigInt d0 ? d0.Mod : ellCoefs.Disc.P;
        if (discP != 0 && discP <= 3)
            throw new GroupException(GroupExceptionType.GroupDef);

        var longCoefs = ec.Flat().ToLongWeierstrassForm();
        var (_, A2l, _, A4l, A6l) = longCoefs.Model;
        var (rl, sl, tl, ul) = longCoefs.TransCoef;
        LongForm = (A2l, A4l, A6l, rl, sl, tl, ul);

        var shortCoefs = ec.Flat().ToShortWeierstrassForm();
        var (_, _, _, A4s, A6s) = shortCoefs.Model;
        var (rs, ss, ts, us) = shortCoefs.TransCoef;
        ShortForm = (A4s, A6s, rs, ss, ts, us);

        Field = $"F{discP}[X,Y]";
        Hash = (Coefs, ShortForm).GetHashCode();
        EllCoefs = ellCoefs;
    }

    public EllGroupSymb(K a1, K a2, K a3, K a4, K a6, EllPoly<K> divPol)
        : this(new EllCoefs<K>(a1, a2, a3, a4, a6), divPol)
    {
    }

    public EllGroupSymb(K a, K b, EllPoly<K> divPol) : this(a.Zero, a.Zero, a.Zero, a, b, divPol)
    {
    }

    public EllGroupSymb(K a, K b, K c, EllPoly<K> divPol) : this(a.Zero, a, a.Zero, b, c, divPol)
    {
    }

    public IEnumerator<EllPt<EllFracPoly<K>>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EllPt<EllFracPoly<K>>>? other) => other?.Hash == Hash;

    public EllPt<EllFracPoly<K>> ConvertToShort(EllPt<EllFracPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return new(x, y);
    }

    public EllPt<EllFracPoly<K>> ConvertFromShort(EllPt<EllFracPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return new(x, y);
    }

    public EllPt<EllFracPoly<K>> ConvertToLong(EllPt<EllFracPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return new(x, y);
    }

    public EllPt<EllFracPoly<K>> ConvertFromLong(EllPt<EllFracPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return new(x, y);
    }

    public int Hash { get; }

    public string Eq => EllCoefs.Eq;

    public string EqShort => $"Elliptic curve short form {EllCoefs.ToShortWeierstrassForm().Eq}";

    public string EqLong => $"Elliptic curve long form {EllCoefs.ToLongWeierstrassForm().Eq}";

    public string Name => a1.IsZero() && a2.IsZero() && a3.IsZero()
        ? $"Ell[{a4},{a6}]({Field})".Replace(" ", "")
        : $"Ell[{a1},{a2},{a3},{a4},{a6}]({Field})".Replace(" ", "");

    public string Field { get; set; }
    public EllFracPoly<K> X { get; }
    public EllFracPoly<K> Y { get; }
    public bool CheckValidity { get; set; } = true;
    public (EllFracPoly<K> a1, EllFracPoly<K> a2, EllFracPoly<K> a3, EllFracPoly<K> a4, EllFracPoly<K> a6) Coefs { get; }

    public EllFracPoly<K> Disc { get; }

    public EllPt<EllFracPoly<K>> Pt
    {
        get
        {
            var lhs = Y.Pow(2) + a1 * X * Y + a3 * Y;
            var rhs = X.Pow(3) + a2 * X.Pow(2) + a4 * X + a6;
            if (rhs.IsZero())
                return new(X, X.Zero);

            return new(X, Y);
        }
    }

    public EllPoly<K> DivPolynomial { get; }

    public (EllFracPoly<K> A, EllFracPoly<K> B, EllFracPoly<K> r, EllFracPoly<K> s, EllFracPoly<K> t, EllFracPoly<K> u)
        ShortForm { get; }

    public (EllFracPoly<K> A, EllFracPoly<K> B, EllFracPoly<K> C, EllFracPoly<K> r, EllFracPoly<K> s, EllFracPoly<K> t,
        EllFracPoly<K> u) LongForm { get; }

    public string ShortFormStr => $"Ell[{ShortForm.A},{ShortForm.B}]({Field})".Replace(" ", "");
    public string LongFormStr => $"Ell[0,{LongForm.A},0,{LongForm.B},{LongForm.C}]({Field})".Replace(" ", "");
    public EllFracPoly<K> a1 => Coefs.a1;
    public EllFracPoly<K> a2 => Coefs.a2;
    public EllFracPoly<K> a3 => Coefs.a3;
    public EllFracPoly<K> a4 => Coefs.a4;
    public EllFracPoly<K> a6 => Coefs.a6;

    public EllFracPoly<K>[] ArrCoefs => [a1, a2, a3, a4, a6];
    public EllCoefs<K> EllCoefs { get; }

    public EllPt<EllFracPoly<K>> this[params ValueType[] us]
    {
        get
        {
            if (us.Length == 2 && us[0] is int x0 && us[1] is int y0)
            {
                var (x1, y1) = (ShortForm.A.One * x0, ShortForm.A.One * y0);
                if (!Contains(x1, y1))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<EllFracPoly<K>>(x1, y1);
            }

            if (us.Length == 2 && us[0] is EllFracPoly<K> x && us[1] is EllFracPoly<K> y)
            {
                if (!Contains(x, y))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<EllFracPoly<K>>(x, y);
            }

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EllPt<EllFracPoly<K>>> GetElements()
    {
        yield return new EllPt<EllFracPoly<K>>();
    }

    public IEnumerable<EllPt<EllFracPoly<K>>> GetGenerators()
    {
        yield return new EllPt<EllFracPoly<K>>();
    }

    public bool ContainsPt(EllPt<EllFracPoly<K>> pt)
    {
        if (pt.IsO)
            return true;

        return Contains(pt.X, pt.Y);
    }

    public bool Contains(EllFracPoly<K> X0, EllFracPoly<K> Y0)
    {
        if (!CheckValidity)
            return true;

        var lhs = Y0 * Y0 + a1 * X0 * Y0 + a3 * Y0;
        var rhs = X0.Pow(3) + a2 * X0 * X0 + a4 * X0 + a6;
        var diff = (lhs - rhs);
        var check = diff.IsZero();
        if (!check)
        {
            var pt = new EllPt<EllFracPoly<K>>(X0, Y0);
            Console.WriteLine(new { pt });
            Console.WriteLine(new { X.Reduction });
            Console.WriteLine(new { lhs, rhs, diff });
        }

        return check;
    }

    public EllPt<EllFracPoly<K>> O => new();
    public EllPt<EllFracPoly<K>> Neutral() => O;

    public EllPt<EllFracPoly<K>> Invert(EllPt<EllFracPoly<K>> P)
    {
        if (P.IsO)
            return P;

        if (!Contains(P.X, P.Y))
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(P.X, -a1 * P.X - a3 - P.Y);
    }

    public EllPt<EllFracPoly<K>> Op(EllPt<EllFracPoly<K>> e1, EllPt<EllFracPoly<K>> e2)
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
        if (!(x1 - x2).IsDivZero())
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) + a1 * alpha - a2 - x2 - x1;
            var y3 = -(alpha * (x3 - x1) + y1) - a1 * x3 - a3;
            return new(x3, y3);
        }
        else
        {
            if (!(y1 + a1 * x2 + a3 + y2).IsDivZero() && !(x1 * a1 + a3 + 2 * y1).IsDivZero())
            {
                var alpha = (3 * x1.Pow(2) - y1 * a1 + 2 * x1 * a2 + a4) / (x1 * a1 + a3 + 2 * y1);
                var x3 = alpha.Pow(2) + a1 * alpha - a2 - 2 * x1;
                var y3 = -(alpha * (x3 - x1) + y1) - a1 * x3 - a3;
                return new(x3, y3);
            }
            else
                return new();
        }
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;

    public EllPt<EllFracPoly<K>> FrobRl(EllPt<EllFracPoly<K>> pt, int n = 1)
    {
        if (pt.IsO)
            return pt;

        var p = pt.X is EllFracPoly<ZnBigInt> x0 ? x0.KOne.Mod : pt.X.P;
        var x = pt.X.FastPow(BigInteger.Pow(p, n));
        var y = pt.Y.FastPow(BigInteger.Pow(p, n));
        
        return new(x, y);
    }
}