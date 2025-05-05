using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup.Integers;

namespace FastGoat.UserGroup.EllCurve;

public struct EllCoefs<K> where K : struct, IFieldElt<K>, IRingElt<K>, IElt<K>
{
    public EllCoefs(K a1, K a2, K a3, K a4, K a6, K r, K s, K t, K u)
    {
        (this.a1, this.a2, this.a3, this.a4, this.a6) = (a1, a2, a3, a4, a6);
        (this.r, this.s, this.t, this.u) = (r, s, t, u);

        b2 = a1 * a1 + 4 * a2;
        b4 = a1 * a3 + 2 * a4;
        b6 = a3 * a3 + 4 * a6;
        b8 = a1 * a1 * a6 - a1 * a3 * a4 + 4 * a2 * a6 + a2 * a3 * a3 - a4 * a4;

        c4 = b2 * b2 - 24 * b4;
        c6 = -b2.Pow(3) + 36 * b2 * b4 - 216 * b6;

        Disc = -b2 * b2 * b8 - 8 * b4.Pow(3) - 27 * b6 * b6 + 9 * b2 * b4 * b6;
        j = c4.Pow(3) / Disc;
    }

    public EllCoefs(K a1, K a2, K a3, K a4, K a6) : this(a1, a2, a3, a4, a6, a1.Zero, a1.Zero, a1.Zero, a1.One)
    {
    }

    public (K a1, K a2, K a3, K a4, K a6) Model => (a1, a2, a3, a4, a6);
    public K[] ArrModel => [a1, a2, a3, a4, a6];
    public (K b2, K b4, K b6, K b8) B_Invariants => (b2, b4, b6, b8);
    public (K c4, K c6) C_Invariants => (c4, c6);
    public K J_Invariant => j;

    public (K a1, K a2, K a3, K a4, K a6, K b2, K b4, K b6, K b8, K c4, K c6, K j) ModelAndInvariants
        => (a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, j);

    public K Disc { get; }
    public K a1 { get; }
    public K a2 { get; }
    public K a3 { get; }
    public K a4 { get; }
    public K a6 { get; }
    public K b2 { get; }
    public K b4 { get; }
    public K b6 { get; }
    public K b8 { get; }
    public K c4 { get; }
    public K c6 { get; }
    public K j { get; }
    public K r { get; }
    public K s { get; }
    public K t { get; }
    public K u { get; }

    public string Eq
    {
        get
        {
            var (x, y) = Ring.Polynomial(Disc, "x", "y").Deconstruct();
            var lhs = y * y + a1 * x * y + a3 * y;
            var rhs = x.Pow(3) + a2 * x * x + a4 * x + a6;
            return $"Elliptic curve {lhs} = {rhs}";
        }
    }

    public string ModelStr => $"[{Model}]".Replace("(", "").Replace(")", "");

    public string B_InvariantsStr => $"B Invariants b2={b2} b4={b4} b6={b6} b8={b8}";

    public string C_InvariantsStr => $"C Invariants c4={c4} c6={c6}";
    public (K r, K s, K t, K u) TransCoef => (r, s, t, u);

    public void Show()
    {
        Console.WriteLine(Eq);
        Console.WriteLine($"Disc = {Disc}");
        Console.WriteLine(B_InvariantsStr);
        Console.WriteLine(C_InvariantsStr);
        Console.WriteLine($"J Invariant j={j}");
        Console.WriteLine($"TransCoef r={r} s={s} t={t} u={u}");
    }

    public EllCoefs<K> Flat() => new(a1, a2, a3, a4, a6);

    public EllPt<K> RevTrans(EllPt<K> P)
    {
        if (P.IsO)
            return P;

        var x = u * u * P.X + r;
        var y = u.Pow(3) * P.Y + s * u * u * P.X + t;
        return new(x, y);
    }

    public EllPt<K> Trans(EllPt<K> P)
    {
        if (P.IsO)
            return P;

        var x = (P.X - r) / (u * u);
        var y = (P.Y - t - s * u * u * x) / u.Pow(3);
        return new(x, y);
    }

    public EllCoefs<K> Transform(K r0, K s0, K t0, K u0)
    {
        var a01 = (a1 + 2 * s0) / u0;
        var a02 = (a2 - s0 * a1 + 3 * r0 - s0 * s0) / u0.Pow(2);
        var a03 = (a3 + r0 * a1 + 2 * t0) / u0.Pow(3);
        var a04 = (a4 - s0 * a3 + 2 * r0 * a2 - (t0 + r0 * s0) * a1 + 3 * r0 * r0 - 2 * s0 * t0) / u0.Pow(4);
        var a06 = (a6 + r0 * a4 + r0 * r0 * a2 + r0.Pow(3) - t0 * a3 - t0 * t0 - r0 * t0 * a1) / u0.Pow(6);

        var r1 = r0 * u.Pow(2) + r;
        var s1 = s0 * u + s;
        var t1 = r0 * s * u.Pow(2) + t0 * u.Pow(3) + t;
        var u1 = u0 * u;

        return new(a01, a02, a03, a04, a06, r1, s1, t1, u1);
    }

    public EllCoefs<K> Transform(int r0, int s0, int t0, int u0)
    {
        var o = Disc.One;
        return Transform(r0 * o, s0 * o, t0 * o, u0 * o);
    }

    public EllCoefs<K> ToLongWeierstrassForm()
    {
        var s0 = -a1 / 2;
        var u0 = s0.IsZero() ? s0.One : s0.One / 2;
        var E = Transform(s0.Zero, s0, s0.Zero, u0);
        var t1 = -E.a3 / 2;
        var u1 = t1.IsZero() ? t1.One : t1.One / 2;
        return E.Transform(t1.Zero, t1.Zero, t1, u1).Simplify();
    }

    public EllCoefs<K> ToShortWeierstrassForm()
    {
        var E = ToLongWeierstrassForm();
        var r1 = -E.a2 / 3;
        var u1 = r1.IsZero() ? r1.One : r1.One / 3;
        return E.Transform(r1, r1.Zero, r1.Zero, u1).Simplify();
    }

    private EllCoefs<K> Simplify()
    {
        if (typeof(K).Equals(typeof(Rational)))
        {
            var pows = new[] { 1, 2, 3, 4, 6 };
            var arr = new[] { a1, a2, a3, a4, a6 }.Select(e => (Rational)(dynamic)e).Select(e => e.Num).ToArray();
            
            var gcd = IntExt.GcdBigInt(arr.Where(e => !e.IsZero).ToArray());
            var gcd2 = IntExt.PrimesDec(gcd)
                .Aggregate(BigInteger.One, (acc, e) => acc * BigInteger.Pow(e.Key, (e.Value - (e.Value % 2)) / 2));
            var dividors = IntExt.DividorsBigInt(gcd2).OrderDescending().ToArray();
            var zip = pows.Zip(arr).ToArray();
            var u0 = dividors.First(d => zip.All(e => e.Second % BigInteger.Pow(d, e.First) == 0));
            var u1 = (K)(dynamic)(new Rational(u0));
            return Transform(u1.Zero, u1.Zero, u1.Zero, u1);
        }
        else
            return this;
    }

    public EllGroup<K> ToEllGroup() => new(this);
}