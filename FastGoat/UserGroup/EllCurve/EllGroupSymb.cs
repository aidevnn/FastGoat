using System.Collections;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;

namespace FastGoat.UserGroup.EllCurve;

public struct EllGroupSymb<K> : IGroup<EllPt<EllPoly<K>>> where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
{
    public EllGroupSymb(EllCoefs<K> ellCoefs, Polynomial<K, Xi> divPol)
    {
        var (y, x) = (Y, X) = GetYX(divPol, ellCoefs);
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

    public EllGroupSymb(K a1, K a2, K a3, K a4, K a6, Polynomial<K, Xi> divPol)
        : this(new EllCoefs<K>(a1, a2, a3, a4, a6), divPol)
    {
    }

    public EllGroupSymb(K a, K b, Polynomial<K, Xi> divPol) : this(a.Zero, a.Zero, a.Zero, a, b, divPol)
    {
    }

    public EllGroupSymb(K a, K b, K c, Polynomial<K, Xi> divPol) : this(a.Zero, a, a.Zero, b, c, divPol)
    {
    }

    public IEnumerator<EllPt<EllPoly<K>>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EllPt<EllPoly<K>>>? other) => other?.Hash == Hash;

    public EllPt<EllPoly<K>> ConvertToShort(EllPt<EllPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return DummyTrackPt(new(x, y));
    }

    public EllPt<EllPoly<K>> ConvertFromShort(EllPt<EllPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, r, s, t, u) = ShortForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return DummyTrackPt(new(x, y));
    }

    public EllPt<EllPoly<K>> ConvertToLong(EllPt<EllPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return DummyTrackPt(new(x, y));
    }

    public EllPt<EllPoly<K>> ConvertFromLong(EllPt<EllPoly<K>> pt)
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
    public EllPoly<K> X { get; }
    public EllPoly<K> Y { get; }
    public bool CheckValidity { get; set; } = true;
    public (EllPoly<K> a1, EllPoly<K> a2, EllPoly<K> a3, EllPoly<K> a4, EllPoly<K> a6) Coefs { get; }

    public EllPoly<K> Disc { get; }

    public EllPt<EllPoly<K>> Pt
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

    public Polynomial<K, Xi> DivPolynomial { get; }

    public (EllPoly<K> A, EllPoly<K> B, EllPoly<K> r, EllPoly<K> s, EllPoly<K> t, EllPoly<K> u)
        ShortForm { get; }

    public (EllPoly<K> A, EllPoly<K> B, EllPoly<K> C, EllPoly<K> r, EllPoly<K> s, EllPoly<K> t,
        EllPoly<K> u) LongForm { get; }

    public string ShortFormStr => $"Ell[{ShortForm.A},{ShortForm.B}]({Field})".Replace(" ", "");
    public string LongFormStr => $"Ell[0,{LongForm.A},0,{LongForm.B},{LongForm.C}]({Field})".Replace(" ", "");
    public EllPoly<K> a1 => Coefs.a1;
    public EllPoly<K> a2 => Coefs.a2;
    public EllPoly<K> a3 => Coefs.a3;
    public EllPoly<K> a4 => Coefs.a4;
    public EllPoly<K> a6 => Coefs.a6;

    public EllPoly<K>[] ArrCoefs => [a1, a2, a3, a4, a6];
    public EllCoefs<K> EllCoefs { get; }

    public EllPt<EllPoly<K>> this[params ValueType[] us]
    {
        get
        {
            if (us.Length == 2 && us[0] is int x0 && us[1] is int y0)
            {
                var (x1, y1) = (ShortForm.A.One * x0, ShortForm.A.One * y0);
                if (!Contains(x1, y1))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<EllPoly<K>>(x1, y1);
            }

            if (us.Length == 2 && us[0] is EllPoly<K> x && us[1] is EllPoly<K> y)
            {
                if (!Contains(x, y))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<EllPoly<K>>(x, y);
            }

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EllPt<EllPoly<K>>> GetElements()
    {
        yield return new EllPt<EllPoly<K>>();
    }

    public IEnumerable<EllPt<EllPoly<K>>> GetGenerators()
    {
        yield return new EllPt<EllPoly<K>>();
    }

    public bool ContainsPt(EllPt<EllPoly<K>> pt)
    {
        if (pt.IsO)
            return true;

        return Contains(pt.X, pt.Y);
    }

    public bool Contains(EllPoly<K> X0, EllPoly<K> Y0)
    {
        if (!CheckValidity)
            return true;

        var lhs = Y0 * Y0 + a1 * X0 * Y0 + a3 * Y0;
        var rhs = X0.Pow(3) + a2 * X0 * X0 + a4 * X0 + a6;
        var diff = (lhs - rhs);
        var check = diff.IsInfty || diff.IsIndeterminate || diff.IsZero();
        if (!check)
        {
            var pt = new EllPt<EllPoly<K>>(X0, Y0);
            Console.WriteLine(new { pt });
            Console.WriteLine(new { X.Reduction });
            Console.WriteLine(new { lhs, rhs, diff });
        }

        return check;
    }

    public EllPt<EllPoly<K>> O => new();
    public EllPt<EllPoly<K>> Neutral() => O;

    public EllPt<EllPoly<K>> Invert(EllPt<EllPoly<K>> P)
    {
        if (P.IsO)
            return P;

        if (!Contains(P.X, P.Y))
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(P.X, -a1 * P.X - a3 - P.Y);
    }

    public EllPt<EllPoly<K>> DummyTrackPt(EllPt<EllPoly<K>> pt)
    {
        if (pt.IsO)
            return pt;

        if (!pt.X.IsNumber || !pt.Y.IsNumber)
        {
            Console.WriteLine("??????????? ######################## ");
            return new();
        }

        if (!Contains(pt.X, pt.Y))
            Console.WriteLine("######################## ??????????? ");

        return pt;
    }

    public EllPt<EllPoly<K>> Op(EllPt<EllPoly<K>> e1, EllPt<EllPoly<K>> e2)
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

    EllPoly<K> FastPowRl(EllPoly<K> a, BigInteger k)
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

    public EllPt<EllPoly<K>> FrobRl(EllPt<EllPoly<K>> pt, int n = 1)
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

    public static EllGroupSymb<K> FromEllGroup(EllGroup<K> E, Polynomial<K, Xi> divPoly)
    {
        return new(E.ToEllCoefs(), divPoly);
    }

    public static EllGroupSymb<K> FromEllGroup(EllGroup<K> E)
    {
        var (_, x) = Ring.Polynomial(E.a1, MonomOrder.Lex, "Y", "X").Deconstruct();
        return new(E.ToEllCoefs(), x.Zero);
    }

    public static EllGroupSymb<K> FromEllCoefs(EllCoefs<K> E, Polynomial<K, Xi> divPoly)
    {
        return new(E, divPoly);
    }

    public static EllGroupSymb<K> FromEllCoefs(EllCoefs<K> E)
    {
        var (_, x) = Ring.Polynomial(E.a1, MonomOrder.Lex, "Y", "X").Deconstruct();
        return new(E, x.Zero);
    }

    public static (EllPoly<K> Y, EllPoly<K> X) GetYX(Polynomial<K, Xi> divPol, EllCoefs<K> ellCoefs)
    {
        var (y, x) = divPol.AllVariables.Deconstruct();
        var (eq, sd) = ellCoefs.ToEllGroup().GetPolynomials(x, y);
        return EllPoly<K>.GetYX(eq, sd, divPol);
    }
}