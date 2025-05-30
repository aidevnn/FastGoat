using System.Collections;
using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.EllCurve;

public struct EllGroupSymb : IGroup<EllPt<EllPoly<ZnInt>>>
{
    public EllGroupSymb(EllCoefs<ZnInt> ellCoefs, KPoly<ZnInt> divPol)
    {
        var (x, y) = (X, Y) = GetXY(divPol, ellCoefs);
        DivPolynomial = divPol;
        var ec = ellCoefs.ToEllCoefs(x);
        Coefs = ec.Model;
        Disc = ec.Disc;

        if (Disc.P != 0 && Disc.P <= 3)
            throw new GroupException(GroupExceptionType.GroupDef);

        var longCoefs = ec.Flat().ToLongWeierstrassForm();
        var (_, A2l, _, A4l, A6l) = longCoefs.Model;
        var (rl, sl, tl, ul) = longCoefs.TransCoef;
        LongForm = (A2l, A4l, A6l, rl, sl, tl, ul);

        var shortCoefs = ec.Flat().ToShortWeierstrassForm();
        var (_, _, _, A4s, A6s) = shortCoefs.Model;
        var (rs, ss, ts, us) = shortCoefs.TransCoef;
        ShortForm = (A4s, A6s, rs, ss, ts, us);

        Field = $"F{Disc.P}[X,Y]";
        Hash = (Coefs, ShortForm).GetHashCode();
        EllCoefs = ellCoefs;
    }

    public EllGroupSymb(ZnInt a1, ZnInt a2, ZnInt a3, ZnInt a4, ZnInt a6, KPoly<ZnInt> divPol)
        : this(new EllCoefs<ZnInt>(a1, a2, a3, a4, a6), divPol)
    {
    }

    public EllGroupSymb(ZnInt a, ZnInt b, KPoly<ZnInt> divPol) : this(a.Zero, a.Zero, a.Zero, a, b, divPol)
    {
    }

    public EllGroupSymb(ZnInt a, ZnInt b, ZnInt c, KPoly<ZnInt> divPol) : this(a.Zero, a, a.Zero, b, c, divPol)
    {
    }

    public IEnumerator<EllPt<EllPoly<ZnInt>>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EllPt<EllPoly<ZnInt>>>? other) => other?.Hash == Hash;

    public EllPt<EllPoly<ZnInt>> ConvertToShort(EllPt<EllPoly<ZnInt>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return DummyTrackPt(new(x, y));
    }

    public EllPt<EllPoly<ZnInt>> ConvertFromShort(EllPt<EllPoly<ZnInt>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return DummyTrackPt(new(x, y));
    }

    public EllPt<EllPoly<ZnInt>> ConvertToLong(EllPt<EllPoly<ZnInt>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return DummyTrackPt(new(x, y));
    }

    public EllPt<EllPoly<ZnInt>> ConvertFromLong(EllPt<EllPoly<ZnInt>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return DummyTrackPt(new(x, y));
    }

    public int Hash { get; }

    public string Eq => EllCoefs.Eq;

    public string EqShort => $"Elliptic curve short form {EllCoefs.ToShortWeierstrassForm().Eq}";

    public string EqLong => $"Elliptic curve long form {EllCoefs.ToLongWeierstrassForm().Eq}";

    public string Name => a1.IsZero() && a2.IsZero() && a3.IsZero()
        ? $"Ell[{a4},{a6}]({Field})".Replace(" ", "")
        : $"Ell[{a1},{a2},{a3},{a4},{a6}]({Field})".Replace(" ", "");

    public string Field { get; set; }
    public EllPoly<ZnInt> X { get; }
    public EllPoly<ZnInt> Y { get; }

    public (EllPoly<ZnInt> a1, EllPoly<ZnInt> a2, EllPoly<ZnInt> a3, EllPoly<ZnInt> a4, EllPoly<ZnInt> a6) Coefs
    {
        get;
    }

    public EllPoly<ZnInt> Disc { get; }

    public EllPt<EllPoly<ZnInt>> Pt
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

    public KPoly<ZnInt> DivPolynomial { get; }

    public (EllPoly<ZnInt> A, EllPoly<ZnInt> B, EllPoly<ZnInt> r, EllPoly<ZnInt> s, EllPoly<ZnInt> t, EllPoly<ZnInt> u)
        ShortForm { get; }

    public (EllPoly<ZnInt> A, EllPoly<ZnInt> B, EllPoly<ZnInt> C, EllPoly<ZnInt> r, EllPoly<ZnInt> s, EllPoly<ZnInt> t,
        EllPoly<ZnInt> u) LongForm { get; }

    public string ShortFormStr => $"Ell[{ShortForm.A},{ShortForm.B}]({Field})".Replace(" ", "");
    public string LongFormStr => $"Ell[0,{LongForm.A},0,{LongForm.B},{LongForm.C}]({Field})".Replace(" ", "");
    public EllPoly<ZnInt> a1 => Coefs.a1;
    public EllPoly<ZnInt> a2 => Coefs.a2;
    public EllPoly<ZnInt> a3 => Coefs.a3;
    public EllPoly<ZnInt> a4 => Coefs.a4;
    public EllPoly<ZnInt> a6 => Coefs.a6;

    public EllPoly<ZnInt>[] ArrCoefs => [a1, a2, a3, a4, a6];
    public EllCoefs<ZnInt> EllCoefs { get; }

    public EllPt<EllPoly<ZnInt>> this[params ValueType[] us]
    {
        get
        {
            if (us.Length == 2 && us[0] is int x0 && us[1] is int y0)
            {
                var (x1, y1) = (ShortForm.A.One * x0, ShortForm.A.One * y0);
                if (!Contains(x1, y1))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<EllPoly<ZnInt>>(x1, y1);
            }

            if (us.Length == 2 && us[0] is EllPoly<ZnInt> x && us[1] is EllPoly<ZnInt> y)
            {
                if (!Contains(x, y))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<EllPoly<ZnInt>>(x, y);
            }

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EllPt<EllPoly<ZnInt>>> GetElements()
    {
        yield return new EllPt<EllPoly<ZnInt>>();
    }

    public IEnumerable<EllPt<EllPoly<ZnInt>>> GetGenerators()
    {
        yield return new EllPt<EllPoly<ZnInt>>();
    }

    public bool ContainsPt(EllPt<EllPoly<ZnInt>> pt)
    {
        if (pt.IsO)
            return true;

        return Contains(pt.X, pt.Y);
    }

    public bool Contains(EllPoly<ZnInt> X0, EllPoly<ZnInt> Y0)
    {
        var lhs = Y0 * Y0 + a1 * X0 * Y0 + a3 * Y0;
        var rhs = X0.Pow(3) + a2 * X0 * X0 + a4 * X0 + a6;
        var diff = (lhs - rhs);
        var check = diff.IsInfty || diff.IsIndeterminate || diff.IsZero();
        if (!check)
        {
            var pt = new EllPt<EllPoly<ZnInt>>(X0, Y0);
            Console.WriteLine(new { pt });
            Console.WriteLine(new { X.Reduction });
            Console.WriteLine(new { lhs, rhs, diff });
        }

        return check;
    }

    public EllPt<EllPoly<ZnInt>> O => new();
    public EllPt<EllPoly<ZnInt>> Neutral() => O;

    public EllPt<EllPoly<ZnInt>> Invert(EllPt<EllPoly<ZnInt>> P)
    {
        if (P.IsO)
            return P;

        if (!Contains(P.X, P.Y))
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(P.X, -a1 * P.X - a3 - P.Y);
    }

    public EllPt<EllPoly<ZnInt>> DummyTrackPt(EllPt<EllPoly<ZnInt>> pt)
    {
        if (pt.IsO)
            return pt;

        if (!pt.X.IsDeterminate || !pt.Y.IsDeterminate)
        {
            Console.WriteLine("??????????? ######################## ");
            return new();
        }

        if (!Contains(pt.X, pt.Y))
            Console.WriteLine("######################## ??????????? ");

        return pt;
    }

    public EllPt<EllPoly<ZnInt>> Op(EllPt<EllPoly<ZnInt>> e1, EllPt<EllPoly<ZnInt>> e2)
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
            return DummyTrackPt(new(x3, y3));
        }
        else
        {
            if (!(y1 + a1 * x2 + a3 + y2).IsDivZero())
            {
                var alpha = (3 * x1.Pow(2) - y1 * a1 + 2 * x1 * a2 + a4) / (x1 * a1 + a3 + 2 * y1);
                var x3 = alpha.Pow(2) + a1 * alpha - a2 - 2 * x1;
                var y3 = -(alpha * (x3 - x1) + y1) - a1 * x3 - a3;
                return DummyTrackPt(new(x3, y3));
            }
            else
                return new();
        }
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;

    EllPoly<ZnInt> FastPowRl(EllPoly<ZnInt> a, BigInteger k)
    {
        if (k == 0)
            return a.One;

        if (k < 0)
            return FastPowRl(a.Inv(), -k);

        var (r, a0, e0) = (a.One, a, k);
        while (e0 > 0)
        {
            if (e0 % 2 == 1)
                r *= a0;

            e0 >>= 1;
            a0 *= a0;
        }

        return r;
    }

    public EllPt<EllPoly<ZnInt>> FrobRl(EllPt<EllPoly<ZnInt>> pt, int n = 1)
    {
        if (pt.IsO)
            return pt;

        var p = pt.X.P;
        var x = FastPowRl(pt.X, BigInteger.Pow(p, n));
        var y = FastPowRl(pt.Y, BigInteger.Pow(p, n));

        if (x.IsIndeterminate || y.IsIndeterminate)
            throw new("Indeterminate exception");
        if (x.IsInfty || y.IsInfty)
            return new();

        return new(x, y);
    }

    public static EllGroupSymb FromEllGroup(EllGroup<ZnInt> E, KPoly<ZnInt> divPoly)
    {
        return new(E.ToEllCoefs(), divPoly);
    }

    public static EllGroupSymb FromEllGroup(EllGroup<ZnInt> E)
    {
        return new(E.ToEllCoefs(), FG.ZPoly(E.a1.P).One);
    }

    public static EllGroupSymb FromEllCoefs(EllCoefs<ZnInt> E)
    {
        return new(E, FG.ZPoly(E.a1.P).One);
    }

    public static EllGroupSymb FromEllCoefs(EllCoefs<ZnInt> E, KPoly<ZnInt> divPoly)
    {
        return new(E, divPoly);
    }

    public static (EllPoly<ZnInt> X, EllPoly<ZnInt> Y) GetXY(KPoly<ZnInt> divPol, EllCoefs<ZnInt> ellCoefs)
    {
        var (x, y, eq, sd) = ellCoefs.ToEllGroup().GetPolynomials();
        var dvp = divPol.Substitute(x);
        return EllPoly<ZnInt>.GetXY(eq, sd, dvp);
    }
}