using System.Collections;
using System.Numerics;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Integers;
using GFelt = FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>;
using RlElt = FastGoat.UserGroup.EllCurve.Frac<FastGoat.UserGroup.EllCurve.Frac<FastGoat.UserGroup.Integers.ZnInt>>;

namespace FastGoat.UserGroup.EllCurve;

public struct EllGroupSymb : IGroup<EllPt<RlElt>>
{
    public EllGroupSymb(EllCoefs<RlElt> ellCoefs, KPoly<ZnInt> divPol)
    {
        Coefs = ellCoefs.Model;
        Disc = ellCoefs.Disc;

        if (!Coefs.a1.IsZero() || !Coefs.a3.IsZero())
            throw new GroupException(GroupExceptionType.GroupDef);

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

        Field = $"F{Disc.P}[X,Y]";

        Hash = (Coefs, ShortForm).GetHashCode();
        (X, Y) = GetXY(Disc.P, 'X', 'Y');
        var dvp = DivPolynomial = divPol;
        EllEqLhs = Y.Pow(2) + a1 * X * Y + a3 * Y;
        EllEqRhs = SimplifyDivPol(X.Pow(3) + a2 * X.Pow(2) + a4 * X + a6, dvp);
        var eq = EllEq = EllEqLhs - EllEqRhs;
        SimplifyPt = pt => SimplifyEllPt(pt, dvp, eq.Num);
    }

    public EllGroupSymb(RlElt a1, RlElt a2, RlElt a3, RlElt a4, RlElt a6, KPoly<ZnInt> divPol)
        : this(new EllCoefs<RlElt>(a1, a2, a3, a4, a6), divPol)
    {
    }

    public EllGroupSymb(RlElt a, RlElt b, KPoly<ZnInt> divPol) : this(a.Zero, a.Zero, a.Zero, a, b, divPol)
    {
    }

    public EllGroupSymb(RlElt a, RlElt b, RlElt c, KPoly<ZnInt> divPol) : this(a.Zero, a, a.Zero, b, c, divPol)
    {
    }

    public IEnumerator<EllPt<RlElt>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<EllPt<RlElt>>? other) => other?.Hash == Hash;

    public EllPt<RlElt> ConvertToShort(EllPt<RlElt> pt)
    {
        if (pt.IsO)
            return pt;

        try
        {
            var (_, _, r, s, t, u) = ShortForm;
            var x = (pt.X - r) / (u * u);
            var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
            return SimplifyPt(new(x, y));
        }
        catch (Exception e)
        {
            var (_, _, r, s, t, u) = ShortForm;
            var x = (pt.X - r) / (u * u);
            var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
            Console.WriteLine(e.Message);
            Console.WriteLine(new { ShortForm });
            Console.WriteLine(new{x});
            Console.WriteLine(new{y});
            return O;
        }
    }

    public EllPt<RlElt> ConvertFromShort(EllPt<RlElt> pt)
    {
        if (pt.IsO)
            return pt;

        try
        {
            var (_, _, r, s, t, u) = ShortForm;
            var x = u * u * pt.X + r;
            var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
            return SimplifyPt(new(x, y));
        }
        catch (Exception e)
        {
            var (_, _, r, s, t, u) = ShortForm;
            var x = u * u * pt.X + r;
            var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
            Console.WriteLine(e.Message);
            Console.WriteLine(new { ShortForm });
            Console.WriteLine(new { x = SimplifyEll(x, DivPolynomial, EllEq.Num) });
            Console.WriteLine(new { y = SimplifyEll(y, DivPolynomial, EllEq.Num) });
            return O;
        }
    }

    public EllPt<RlElt> ConvertToLong(EllPt<RlElt> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = (pt.X - r) / (u * u);
        var y = (pt.Y - t - s * u * u * x) / u.Pow(3);
        return new(x, y);
    }

    public EllPt<RlElt> ConvertFromLong(EllPt<RlElt> pt)
    {
        if (pt.IsO)
            return pt;

        var (_, _, _, r, s, t, u) = LongForm;
        var x = u * u * pt.X + r;
        var y = u.Pow(3) * pt.Y + s * u * u * pt.X + t;
        return new(x, y);
    }

    public Func<EllPt<RlElt>, EllPt<RlElt>> SimplifyPt { get; }

    public int Hash { get; }

    public string Eq => $"Elliptic curve {EllEqLhs} = {EllEqRhs}";

    public string EqShort => $"Elliptic curve short form {Y * Y} = {X.Pow(3) + ShortForm.A * X + ShortForm.B}";

    public string EqLong =>
        $"Elliptic curve long form {Y * Y} = {X.Pow(3) + LongForm.A * X * X + LongForm.B * X + LongForm.C}";

    public string Name => a1.IsZero() && a2.IsZero() && a3.IsZero()
        ? $"Ell[{a4},{a6}]({Field})".Replace(" ", "")
        : $"Ell[{a1},{a2},{a3},{a4},{a6}]({Field})".Replace(" ", "");

    public string Field { get; set; }
    public RlElt X { get; }
    public RlElt Y { get; }
    public (RlElt a1, RlElt a2, RlElt a3, RlElt a4, RlElt a6) Coefs { get; }
    public RlElt Disc { get; }
    public RlElt EllEqLhs { get; }
    public RlElt EllEqRhs { get; }
    public RlElt EllEq { get; }

    public EllPt<RlElt> Pt
    {
        get
        {
            if (EllEqRhs.IsZero())
                return new(X, Y.Zero);

            return new(X, Y);
        }
    }

    public KPoly<ZnInt> DivPolynomial { get; }
    public (RlElt A, RlElt B, RlElt r, RlElt s, RlElt t, RlElt u) ShortForm { get; }
    public (RlElt A, RlElt B, RlElt C, RlElt r, RlElt s, RlElt t, RlElt u) LongForm { get; }
    public string ShortFormStr => $"Ell[{ShortForm.A},{ShortForm.B}]({Field})".Replace(" ", "");
    public string LongFormStr => $"Ell[0,{LongForm.A},0,{LongForm.B},{LongForm.C}]({Field})".Replace(" ", "");
    public RlElt a1 => Coefs.a1;
    public RlElt a2 => Coefs.a2;
    public RlElt a3 => Coefs.a3;
    public RlElt a4 => Coefs.a4;
    public RlElt a6 => Coefs.a6;

    public RlElt[] ArrCoefs => [a1, a2, a3, a4, a6];

    public EllPt<RlElt> this[params ValueType[] us]
    {
        get
        {
            if (us.Length == 2 && us[0] is int x0 && us[1] is int y0)
            {
                var (x1, y1) = (ShortForm.A.One * x0, ShortForm.A.One * y0);
                if (!Contains(x1, y1))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<RlElt>(x1, y1);
            }

            if (us.Length == 2 && us[0] is RlElt x && us[1] is RlElt y)
            {
                if (!Contains(x, y))
                    throw new GroupException(GroupExceptionType.GroupDef);

                return new EllPt<RlElt>(x, y);
            }

            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

    public IEnumerable<EllPt<RlElt>> GetElements()
    {
        yield return new EllPt<RlElt>();
    }

    public IEnumerable<EllPt<RlElt>> GetGenerators()
    {
        yield return new EllPt<RlElt>();
    }

    public bool ContainsPt(EllPt<RlElt> pt)
    {
        if (pt.IsO)
            return true;

        return Contains(pt.X, pt.Y);
    }

    public bool Contains(RlElt X0, RlElt Y0)
    {
        var lhs = Y0 * Y0 + a1 * X0 * Y0 + a3 * Y0;
        var rhs = X0.Pow(3) + a2 * X0 * X0 + a4 * X0 + a6;
        var diff = SimplifyEll(lhs - rhs, DivPolynomial, EllEq.Num);
        return diff.IsInfinity || diff.IsZero();
    }

    public EllPt<RlElt> O => new();
    public EllPt<RlElt> Neutral() => O;

    public EllPt<RlElt> Invert(EllPt<RlElt> P)
    {
        if (P.IsO)
            return P;

        if (!Contains(P.X, P.Y))
            throw new GroupException(GroupExceptionType.GroupDef);

        return SimplifyPt(new(P.X, -a1 * P.X - a3 - P.Y));
    }

    public EllPt<RlElt> Op(EllPt<RlElt> e1, EllPt<RlElt> e2)
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
        if (!SimplifyEll(x1 - x2, DivPolynomial, EllEq.Num).IsZero())
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) + a1 * alpha - a2 - x2 - x1;
            var y3 = x3 * alpha - x1 * alpha + y1;
            return SimplifyPt(new(x3, -a1 * x3 - a3 - y3));
        }
        else
        {
            if (!SimplifyEll(y1 + a1 * x2 + a3 + y2, DivPolynomial, EllEq.Num).IsZero())
            {
                var alpha = (3 * x1.Pow(2) - y1 * a1 + 2 * x1 * a2 + a4) / (x1 * a1 + a3 + 2 * y1);
                var x3 = alpha.Pow(2) + a1 * alpha - a2 - 2 * x1;
                var y3 = alpha * (x3 - x1) + y1;
                return SimplifyPt(new(x3, -a1 * x3 - a3 - y3));
            }
            else
                return new();
        }
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;

    RlElt FastPowRl(RlElt a, BigInteger k)
    {
        if (k == 0)
            return a.One;

        if (k < 0)
            return FastPowRl(a.Inv(), -k);

        var (r, a0, e0) = (a.One, a, k);
        while (e0 > 0)
        {
            if (e0 % 2 == 1)
                r = SimplifyEll(r * a0, DivPolynomial, EllEq.Num);

            e0 >>= 1;
            a0 = SimplifyEll(a0 * a0, DivPolynomial, EllEq.Num);
        }

        return r;
    }

    public EllPt<RlElt> FrobRl(EllPt<RlElt> pt, int n = 1)
    {
        if (pt.IsO)
            return pt;

        var p = pt.X.P;
        var x = FastPowRl(pt.X, BigInteger.Pow(p, n));
        var y = FastPowRl(pt.Y, BigInteger.Pow(p, n));

        if (x.IsIndeterminate || y.IsIndeterminate)
            throw new("Indeterminate exception");
        if (x.IsInfinity || y.IsInfinity)
            return new();

        return new(x, y);
    }

    (GFelt num, GFelt denom) Substitute(RlElt P1, EllPt<GFelt> P2)
    {
        if (P1.Num.Coefs.Any(c => c.Denom.Substitute(P2.X).IsZero()) ||
            P1.Denom.Coefs.Any(c => c.Denom.Substitute(P2.X).IsZero()))
            return (P2.X.Zero, P2.X.Zero);

        var num = P1.Num.Coefs.Select((c, i) => P2.Y.Pow(i) * c.Num.Substitute(P2.X) / c.Denom.Substitute(P2.X))
            .ToVec().Sum();
        var denom = P1.Denom.Coefs.Select((c, i) => P2.Y.Pow(i) * c.Num.Substitute(P2.X) / c.Denom.Substitute(P2.X))
            .ToVec().Sum();

        return (num, denom);
    }

    public EllPt<GFelt> SubstitutePt(EllPt<RlElt> P1, EllPt<GFelt> P2)
    {
        if (P1.IsO)
            return new();

        var x = Substitute(P1.X, P2);
        var y = Substitute(P1.Y, P2);
        if (x.denom.IsZero() || y.denom.IsZero())
            return new();

        return new(x.num / x.denom, y.num / y.denom);
    }

    public static KPoly<ZnInt> SimplifyDivPol(KPoly<ZnInt> P, KPoly<ZnInt> divPoly)
    {
        var (quo, rem) = P.Div(divPoly);
        if (rem.IsZero())
            return rem;

        if (quo.IsZero())
        {
            var gcd = Ring.FastGCD(divPoly, P);
            if (gcd.Degree > 1)
                return P.Zero;

            return P;
        }

        return SimplifyDivPol(rem, divPoly);
    }

    public static Frac<ZnInt> SimplifyDivPol(Frac<ZnInt> P, KPoly<ZnInt> divPoly)
    {
        var num = SimplifyDivPol(P.Num, divPoly);
        var denom = SimplifyDivPol(P.Denom, divPoly);
        return new(num, denom);
    }

    public static RlElt SimplifyDivPol(RlElt P, KPoly<ZnInt> divPoly)
    {
        var num = new KPoly<Frac<ZnInt>>(P.x, P.Num.KZero,
            P.Num.Coefs.Select(n0 => SimplifyDivPol(n0, divPoly)).ToArray());
        var denom = new KPoly<Frac<ZnInt>>(P.x, P.Denom.KZero,
            P.Denom.Coefs.Select(n0 => SimplifyDivPol(n0, divPoly)).ToArray());
        return new(num, denom);
    }

    public static KPoly<Frac<ZnInt>> SimplifyEll(KPoly<Frac<ZnInt>> P, KPoly<Frac<ZnInt>> EllEq)
    {
        var (quo, rem) = P.Div(EllEq);
        if (quo.IsZero())
            return P;

        return SimplifyEll(rem, EllEq);
    }

    public static RlElt SimplifyEll(RlElt P, KPoly<ZnInt> divPoly, KPoly<Frac<ZnInt>> EllEq)
    {
        var num = SimplifyEll(P.Num, EllEq);
        var denom = SimplifyEll(P.Denom, EllEq);
        var P2 = new RlElt(num, denom);
        var P3 = SimplifyDivPol(P2, divPoly);
        return P3;
    }

    public static EllPt<RlElt> SimplifyEllPt(EllPt<RlElt> Pt, KPoly<ZnInt> divPoly, KPoly<Frac<ZnInt>> EllEq)
    {
        var x = SimplifyEll(Pt.X, divPoly, EllEq);
        var y = SimplifyEll(Pt.Y, divPoly, EllEq);
        if (x.IsIndeterminate || y.IsIndeterminate)
            throw new("Indeterminate exception");
        if (x.IsInfinity || y.IsInfinity)
            return new();

        return new(x, y);
    }

    public static EllGroupSymb FromEllGroup(EllGroup<ZnInt> E, KPoly<ZnInt> divPoly)
    {
        var Y = GetXY(divPoly.P, 'X', 'Y').Y;
        var (a1, a2, a3, a4, a6) = E.Coefs;
        var (A1, A2, A3, A4, A6) = (a1 * Y.KOne * Y.One, a2 * Y.KOne * Y.One, a3 * Y.KOne * Y.One, a4 * Y.KOne * Y.One,
            a6 * Y.KOne * Y.One);
        return new(A1, A2, A3, A4, A6, divPoly);
    }
    
    public static Frac<Rational> FracQ(char x = 'x') => new(FG.QPoly(x));
    public static Frac<ZnInt> FracZnInt(int p, char x = 'x') => new(FG.ZPoly(p, x));

    public static (Frac<Frac<Rational>> X, Frac<Frac<Rational>> Y) BiVariateFracQ(char x = 'x', char y = 'y')
    {
        var X = FracQ(x);
        var Y = new Frac<Frac<Rational>>(y, X);
        return (X * Y.One, Y);
    }

    public static (Frac<Frac<ZnInt>> X, Frac<Frac<ZnInt>> Y) GetXY(int p, char x = 'x', char y = 'y')
    {
        var X = FracZnInt(p, x);
        var Y = new Frac<Frac<ZnInt>>(y, X);
        return (X * Y.One, Y);
    }

}