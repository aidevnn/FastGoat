using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static ZnInt ToZnInt(this Rational e, int p) =>
        new ZnInt(p, (int)BigInteger.Remainder(e.Num, p)) / new ZnInt(p, (int)BigInteger.Remainder(e.Denom, p));

    public static ZnBInt ToZnBInt(this Rational e, Modulus p) =>
        new ZnBInt(p, BigInteger.Remainder(e.Num, p.Mod)) / new ZnBInt(p, BigInteger.Remainder(e.Denom, p.Mod));

    public static ZnBigInt ToZnBigInt(this Rational e, BigInteger p) =>
        new ZnBigInt(p, BigInteger.Remainder(e.Num, p)) / new ZnBigInt(p, BigInteger.Remainder(e.Denom, p));

    public static EPoly<ZnInt> ToGF(this Rational e, BigInteger q, char a = 'a')
    {
        var x = FqX(q, a);
        return x.One * e.ToZnInt(x.P);
    }

    public static EPoly<ZnInt> ToGF(this EPoly<ZnInt> e, BigInteger q, char a = 'a')
    {
        if (!e.F.IsInConwayPolynomialsDB())
            throw new($"[{e.P},{e.F.Degree}:{e.F.Glue(",")}] is not in local Conway Polynomials Database");
        
        var x = FqX(q, a);
        if (x.P != e.P)
            throw new();
        
        var q0 = BigInteger.Pow(e.P, e.F.Degree);
        if (x.F.Degree % e.F.Degree != 0)
            throw new();

        var d = (q - 1) / (q0 - 1);
        return e.Substitute(x.FastPow(d));
    }

    public static Rational Comb(int k, int n)
    {
        var c = Rational.KOne();
        if (k > n || k < 0)
            return c;

        var num = c;
        var denom = c;
        k = 2 * k > n ? n - k : k;
        for (int i = 0; i < k; i++)
        {
            num *= (n - i);
            denom *= (i + 1);
        }

        return num / denom;
    }

    public static ((int, int)[]num, (int, int)[] denom) Decomp(this Rational r)
    {
        var num = IntExt.PrimesDec(r.Num).Select(e => (e.Key, e.Value)).ToArray();
        var denom = IntExt.PrimesDec(r.Denom).Select(e => (e.Key, e.Value)).ToArray();
        return (num, denom);
    }

    public static KPoly<ZnInt> ZPoly(int p, char x = 'x') => new KPoly<ZnInt>(x, ZnInt.ZnZero(p)).X;
    public static KPoly<ZnBInt> ZbPoly(int p, char x = 'x') => new KPoly<ZnBInt>(x, ZnBInt.ZnZero(p)).X;
    public static KPoly<ZnBInt> ZbPoly(int p, int o, char x = 'x') => new KPoly<ZnBInt>(x, new ZnBInt(new(p, o), 0)).X;
    public static KPoly<ZnBigInt> ZlPoly(int p, char x = 'x') => new KPoly<ZnBigInt>(x, ZnBigInt.ZnZero(p)).X;
    public static KPoly<Rational> QPoly(char x = 'x') => new(x);
    public static KPoly<Cplx> CplxPoly(char x = 'x') => new(x);
    public static KPoly<BigCplx> BCplxPoly(int o = 40, char x = 'x') => new KPoly<BigCplx>(x, BigCplx.BcZero(o)).X;
    public static KPoly<BigReal> BRealPoly(int o = 40, char x = 'x') => new KPoly<BigReal>(x, BigReal.BrZero(o)).X;

    public static Cplx NSolve(KPoly<Cplx> P, double epsilon = 1e-14, int maxLoop = 200)
    {
        var dP = P.Derivative;
        var ai = Cplx.CZero;
        var aj = new Cplx(Double.Pi + Complex.ImaginaryOne);
        var i = 0;
        do
        {
            ai = aj;
            aj = ai - (P.Substitute(ai) / dP.Substitute(ai)); // Newton iteration
        } while ((ai - aj).NormInf > epsilon && i++ < maxLoop);

        return aj;
    }

    public static Cplx[] NRoots(KPoly<Cplx> P, double epsilon = 1e-14, int maxLoop = 200)
    {
        var P0 = P;
        var roots = new List<Cplx>();
        while (P0.Degree > 0)
        {
            var a0 = NSolve(P0, epsilon, maxLoop);
            roots.Add(a0);
            P0 /= P0.X - a0;
        }

        return roots.ToArray();
    }

    public static BigCplx NSolve(KPoly<BigCplx> P, int maxLoop = 200)
    {
        var o = P.KZero.O;
        var aj = new BigCplx(BigReal.Pi(o), BigReal.E(o));
        return NSolve(P, aj, maxLoop);
    }

    public static BigCplx NSolve(KPoly<BigCplx> P, BigCplx a0, int maxLoop = 200)
    {
        var dP = P.Derivative;
        var i = 0;
        BigCplx ai;
        var aj = a0;
        do
        {
            ai = aj;
            var dpai = dP.Substitute(ai);
            if (dpai.IsZero())
                break;
            aj = ai - (P.Substitute(ai) / dpai); // Newton iteration
        } while (!(ai - aj).IsZero() && i++ < maxLoop);

        return aj;
    }

    public static BigCplx[] NRoots(KPoly<BigCplx> P, int maxLoop = 200)
    {
        var O1 = P.KZero.O;
        var O2 = O1 / 2;
        var P0 = P;
        var roots = new List<BigCplx>();
        var (pi, e) = (BigReal.Pi(O2), BigReal.E(O2));
        var a0 = new BigCplx(pi, e);

        while (P0.Degree > 0)
        {
            var P1 = P0.ToBcPoly(O2);
            a0 = NSolve(P1, BigCplx.FromBigReal(a0.ImaginaryPart / 2, a0.RealPart / 2), maxLoop / 2);
            if (!P1.Substitute(a0).ToBigCplx(O2 / 2).IsZero())
            {
                a0 = new BigCplx(pi, e);
                continue;
            }

            var a1 = NSolve(P, a0.ToBigCplx(O1), maxLoop);
            roots.Add(a1);
            P0 /= P0.X - a1;
        }

        return roots.ToArray();
    }

    public static bool IsInConwayPolynomialsDB<K>(this KPoly<K> P) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var (n, p) = (P.Degree, P.P);
        return PolynomExt.AllConwayPolys[p][n].Index().All(e => (P[e.Index] - e.Item).IsZero());
    }

    public static KPoly<ZnInt> ToZnPoly(this KPoly<Rational> P, int p) =>
        new(P.x, ZnInt.ZnZero(p), P.Coefs.Select(c => c.ToZnInt(p)).TrimSeq().ToArray());

    public static KPoly<ZnBigInt> ToZnBigIntPoly(this KPoly<Rational> P, BigInteger p) =>
        new(P.x, ZnBigInt.ZnZero(p), P.Coefs.Select(c => c.ToZnBigInt(p)).TrimSeq().ToArray());

    public static KPoly<EPoly<ZnInt>> ToGF(this KPoly<Rational> P, BigInteger q, char x = 'a') =>
        P.Coefs.Select(c => c.ToGF(q, x)).ToKPoly(P.x);

    public static KPoly<EPoly<ZnInt>> ToGF(this KPoly<EPoly<ZnInt>> P, BigInteger q, char x = 'a') =>
        P.Coefs.Select(c => c.ToGF(q, x)).ToKPoly(P.x);

    public static KPoly<Rational> ToRationalPoly(this KPoly<ZnInt> P) =>
        new(P.x, Rational.KOne(), P.Coefs.Select(c => c.K * Rational.KOne()).TrimSeq().ToArray());

    public static KPoly<Rational> RoundPoly(this KPoly<Rational> P) =>
        new(P.x, Rational.KOne(), P.Coefs.Select(c => c.RoundEven).TrimSeq().ToArray());

    public static KPoly<Rational> TruncPoly(this KPoly<Rational> P) =>
        new(P.x, Rational.KOne(), P.Coefs.Select(c => c.Trunc).TrimSeq().ToArray());

    public static KPoly<Rational> FloorPoly(this KPoly<Rational> P) =>
        new(P.x, Rational.KOne(), P.Coefs.Select(c => c.Floor).TrimSeq().ToArray());

    public static KPoly<Cplx> ToCPoly(this KPoly<Rational> P) =>
        new(P.x, Cplx.CZero, P.Coefs.Select(c => c * Cplx.COne).TrimSeq().ToArray());

    public static KPoly<BigReal> ToBrPoly(this KPoly<Rational> P, int O = 40) =>
        new(P.x, BigReal.BrZero(O), P.Coefs.Select(c => BigReal.BrOne(O) * BigReal.FromRational(c, O)).TrimSeq()
            .ToArray());

    public static KPoly<BigCplx> ToBcPoly(this KPoly<Rational> P, int O = 40) =>
        new(P.x, BigCplx.BcZero(O), P.Coefs.Select(c => BigCplx.BcOne(O) * BigCplx.FromRational(c, O)).TrimSeq()
            .ToArray());

    public static KPoly<BigCplx> ToRoundBcPoly(this KPoly<BigCplx> P, int d = 40) =>
        new(P.x, BigCplx.Round(P.KZero, d), P.Coefs.Select(c => BigCplx.Round(c, d)).TrimSeq().ToArray());

    public static KPoly<BigCplx> ToBcPoly(this KPoly<BigCplx> P, int O = 40) =>
        new(P.x, BigCplx.BcZero(O), P.Coefs.Select(c => c.ToBigCplx(O)).TrimSeq().ToArray());

    public static KPoly<Rational> ToIntPoly(this KPoly<BigCplx> P, int err = 4)
    {
        var one = P.KOne;
        if (err < 4 || one.O < 2 * err)
            throw new($"Digits err {err} must be >=4 and <=O/2");

        Rational Cv(BigCplx r)
        {
            if (!r.ImaginaryPart.ToBigReal(one.O - err).IsZero())
                throw new("Only real numbers allowed");

            if (r.ToBigCplx(one.O - err).IsZero())
                return Rational.KZero();

            return r.RealPart.ToBigReal(one.O - err).ToRational.RoundEven;
        }

        return new(P.x, Rational.KZero(), P.Coefs.Select(Cv).TrimSeq().ToArray());
    }

    public static KPoly<Rational> ToIntPoly(this KPoly<BigReal> P)
    {
        var one = P.KOne;

        Rational Cv(BigReal r)
        {
            if (r.IsZero())
                return Rational.KZero();

            return BigReal.Round(r, one.O - 4).ToRational;
        }

        return new(P.x, Rational.KZero(), P.Coefs.Select(Cv).TrimSeq().ToArray());
    }

    public static KPoly<Cplx> ToCPoly(this KPoly<EPoly<Rational>> P, Cplx e) =>
        new(P.x, Cplx.CZero, P.Coefs.Select(c => c.Poly.Substitute(e)).TrimSeq().ToArray());

    public static KPoly<Rational> ToAbsKPoly(this KPoly<Rational> P) =>
        new(P.x, P.KZero, P.Coefs.Select(c => c.Absolute).TrimSeq().ToArray());

    public static KPoly<Rational> Primitive(this KPoly<Rational> P)
    {
        if (P.Degree == 0)
            return P * new Rational(P[0].Denom);

        var lcm = new Rational(IntExt.LcmBigInt(P.Coefs.Select(e => e.Denom).Distinct().ToArray()));
        var P0 = P * lcm; // removes denominators
        var gcd = new Rational(
            IntExt.GcdBigInt(P0.Coefs.Select(e => BigInteger.Abs(e.Num)).Where(e => !e.IsZero).Distinct().ToArray()));
        return P0 / gcd * P.LT.Sign; // makes primitive
    }

    public static KPoly<Rational> ZPoly(this KPoly<Rational> P)
    {
        if (P.Degree == 0)
            return P * new Rational(P[0].Denom);

        var lcm = new Rational(IntExt.LcmBigInt(P.Coefs.Select(e => e.Denom).Distinct().ToArray()));
        return P * lcm; // removes denominators
    }

    public static KPoly<Rational> PrimitiveZPoly(this KPoly<Rational> P) => P.ZPoly().Primitive();

    public static Cnf Substitute(this KPoly<Rational> P, Cnf c)
    {
        return new(c.N, P.Substitute(c.E));
    }

    public static KPoly<Cnf> ToCnfPoly(this KPoly<EPoly<Rational>> P)
    {
        var Fs = P.Coefs.Select(c => c.F).Distinct().ToArray();
        var keys = CyclotomicPolynomials.Where(e => e.Value.Coefs.SequenceEqual(Fs[0].Coefs)).Take(1).ToArray();
        if (Fs.Length != 1 || keys.Length != 1)
        {
            keys.Println("Keys");
            throw new($"Fs:{Fs.Glue(", ")}");
        }

        var N = keys[0].Key;
        var arr = P.Coefs.Select(c => new Cnf(N, c)).TrimSeq().ToArray();
        return new('X', Cnf.CnfZero, arr);
    }

    public static KPoly<Cnf> ToCnfPoly(this KPoly<EPoly<Rational>> P, int N)
    {
        var pol = CyclotomicPolynomial(N);
        var Fs = P.Coefs.Select(c => c.F).Distinct().ToArray();
        if (Fs.Length != 1 || !Fs[0].Coefs.SequenceEqual(pol.Coefs))
            throw new($"Fs:{Fs.Glue(", ")}");

        var arr = P.Coefs.Select(c => new Cnf(N, c)).TrimSeq().ToArray();
        return new('X', Cnf.CnfZero, arr);
    }

    public static KPoly<Cnf> ToCnfPoly(this KPoly<Cnf> P, int N)
    {
        return new('X', Cnf.CnfZero, P.Coefs.Select(c => c.ToCnfN(N)).TrimSeq().ToArray());
    }

    public static Cnf ToCnfN(this Cnf c, int N)
    {
        var lcm = (N * c.N) / (IntExt.Gcd(N, c.N));
        var a = CyclotomicEPoly(lcm);
        var a1 = a.Pow(lcm / c.N);
        return new(lcm, c.E.Substitute(a1));
    }

    public static KPoly<EPoly<Rational>> ToEPolyX(this KPoly<Cnf> P, char x = 'X', char y = 'y')
    {
        var Ns = P.Coefs.Select(c => c.N).ToArray();
        var lcm = IntExt.Lcm(Ns);
        var a = CyclotomicEPoly(lcm);
        var arr = P.Coefs.Select(c => c.E.Substitute(a.Pow(lcm / c.N))).TrimSeq().ToArray();
        return new(x, a.Zero, arr);
    }

    public static Rational NormB(this KPoly<Rational> P, int b)
    {
        if (b < 1)
            throw new();

        return P.Coefs.Aggregate(P.KZero, (acc, c) => acc + c.Absolute.Pow(b));
    }

    public static Rational NormInf(this KPoly<Rational> P)
    {
        return P.Coefs.Max(c => c.Absolute);
    }

    public static int NbCoeffs<K>(this KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return f.Coefs.Count(c => !c.IsZero());
    }

    public static KPoly<K> KPoly<K>(char x, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new KPoly<K>(x, scalar).X;
    }

    public static FracPoly<K> KFracPoly<K>(char x, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new FracPoly<K>(x, scalar).X;
    }

    public static FracPoly<K> KFracPoly<K>(KPoly<K> p) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new FracPoly<K>(p).X;
    }

    public static FracPoly<FracPoly<K>> ToFrac<K>(this KPoly<KPoly<K>> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var x = KFracPoly(f.KOne.X);
        var y = KFracPoly(f.x, x.One);
        return f.Coefs.Select((cy, i) => cy.Substitute(x) * y.Pow(i)).ToVec().Sum();
    }

    public static T Substitute<K, T>(this FracPoly<FracPoly<K>> f, T x, T y)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
        where T : struct, IElt<T>, IRingElt<T>, IFieldElt<T>, IModuleElt<K, T>, IVsElt<K, T>
    {
        var num = f.Num.Coefs.Select((cy, i) => y.Pow(i) * cy.Substitute(x)).ToVec().Sum();
        var denom = f.Denom.Coefs.Select((cy, i) => y.Pow(i) * cy.Substitute(x)).ToVec().Sum();
        return num / denom;
    }

    public static EPoly<K> EPoly<K>(KPoly<K> f) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(f);
    }

    public static EPoly<K> EPoly<K>(KPoly<K> f, char x) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new(new KPoly<K>(x, f.KZero, f.Coefs));
    }

    public static EPoly<Rational> CyclotomicEPoly(int n, char x = 'a') => FG.EPoly(FG.CyclotomicPolynomial(n), x);

    public static KPoly<K> KPoly<K>(char x, params K[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var coefs0 = coefs.TrimSeq().ToArray();
        if (coefs0.Length < 2)
            throw new GroupException(GroupExceptionType.GroupDef);

        return new KPoly<K>(x, coefs[0].Zero, coefs0);
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(KPoly<K> f, char c, char x = 'X')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var c0 = EPoly(f, c);
        var x0 = KPoly(x, c0);
        return (x0, c0 * x0.One);
    }

    public static (KPoly<EPoly<K>> X, EPoly<K> c) EPolyXc<K>(KPoly<K> f, char c, char x = 'X')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var c0 = EPoly(f, c);
        var x0 = KPoly(x, c0);
        return (x0, c0);
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(KPoly<K> f)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPolyXC(f, f.x);
    }

    public static (KPoly<EPoly<K>> X, KPoly<EPoly<K>> c) EPolyXC<K>(K scalar, char a, char b, params int[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPolyXC(coefs.ToKPoly(scalar, a), a, b);
    }

    public static EPoly<K> EPoly<K>(K scalar, char x, params int[] coefs)
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return EPoly(coefs.ToKPoly(scalar, x));
    }

    public static KPoly<K> ResMod<K>(this KPoly<K> P, KPoly<K> F) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
        => P.Div(F).rem;

    public static KPoly<K> ToKPoly<K>(this IEnumerable<K> coefs, char x = 'x')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        var arr0 = coefs.ToArray();
        if (arr0.Length == 0)
            throw new GroupException(GroupExceptionType.GroupDef);
        var arr = arr0.TrimSeq().ToArray();
        return new(x, arr0[0], arr.ToArray());
    }

    public static KPoly<K> ToKPoly<K>(this IEnumerable<int> coefs, K scalar, char x = 'x')
        where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
        => coefs.Select(i => i * scalar.One).ToKPoly(x);

    public static KPoly<Rational> CoefsMod(this KPoly<Rational> P, Rational Q)
        => P.Coefs.Select(c => c.Mod(Q)).ToKPoly();

    public static KPoly<Rational> CoefsMod(this KPoly<Rational> P, int Q)
        => P.Coefs.Select(c => c.Mod(Q)).ToKPoly();

    public static KPoly<Rational> ClosestModulusTo(this KPoly<Rational> source, KPoly<Rational> destination,
        Rational mod)
    {
        var deg = int.Max(source.Degree, destination.Degree);
        return (deg + 1).SeqLazy().Select(i => source[i].ClosestModulusTo(destination[i], mod)).ToKPoly();
    }

    public static KPoly<Rational> ClosestModulusTo(this KPoly<Rational> source, KPoly<Rational> destination, int mod)
    {
        var deg = int.Max(source.Degree, destination.Degree);
        return (deg + 1).SeqLazy().Select(i => source[i].ClosestModulusTo(destination[i], mod)).ToKPoly();
    }

    public static KPoly<Rational> CoefsModSigned(this KPoly<Rational> P, Rational Q)
        => P.Coefs.Select(c => c.Signed(Q)).ToKPoly();

    public static KPoly<Rational> CoefsModSigned(this KPoly<Rational> P, int Q)
        => P.Coefs.Select(c => c.Signed(Q)).ToKPoly();

    public static KPoly<Rational> ResMod(this KPoly<Rational> P, KPoly<Rational> F, Rational Q)
        => P.Div(F).rem.CoefsMod(Q);

    public static KPoly<Rational> ResMod(this KPoly<Rational> P, KPoly<Rational> F, int Q)
        => P.Div(F).rem.CoefsMod(Q);

    public static KPoly<Rational> ResModSigned(this KPoly<Rational> P, KPoly<Rational> F, Rational Q)
        => P.Div(F).rem.CoefsModSigned(Q);

    public static KPoly<Rational> ResModSigned(this KPoly<Rational> P, KPoly<Rational> F, int Q)
        => P.Div(F).rem.CoefsModSigned(Q);

    public static EPoly<Rational> EQPoly(char x, params int[] coefs) => EPoly(Rational.KZero(), x, coefs);

    public static FracPoly<Rational> QFracPoly(char x = 'x') => KFracPoly(x, Rational.KZero());

    public static FracPoly<ZnInt> ZFracPoly(int p, char x = 'x') => KFracPoly(x, ZnInt.ZnZero(p));

    public static (KPoly<FracPoly<ZnInt>> x, FracPoly<ZnInt> t) FpT_Poly(int p, (char x, char t) xt)
    {
        var t = ZPoly(p, xt.t);
        var t0 = new FracPoly<ZnInt>(t);
        return (KPoly(xt.x, t0), t0);
    }

    public static (KPoly<FracPoly<Rational>> X, FracPoly<Rational> T) QT_Poly(char x, char t)
    {
        var t0 = QPoly(t);
        var t1 = new FracPoly<Rational>(t0);
        return (KPoly(x, t1), t1);
    }

    public static (KPoly<FracPoly<ZnInt>> x, FracPoly<ZnInt> t) FpT_Poly(int p) => FpT_Poly(p, ('x', 't'));

    public static (KPoly<EPoly<ZnInt>> x, EPoly<ZnInt> a) FqX_Poly(int q, (char x, char a) xa)
    {
        var a = FqX(q, xa.a);
        return (KPoly(xa.x, a), a);
    }

    public static (KPoly<EPoly<ZnInt>> x, EPoly<ZnInt> a) FqX_Poly(int q) => FqX_Poly(q, ('x', 'a'));

    public static BigInteger GLnqOrder(int n, int q) => BigInteger.Pow(q, n * (n - 1) / 2) * n.Range()
        .Select(k => BigInteger.Pow(q, n - k) - 1)
        .Aggregate(BigInteger.One, (acc, pk) => acc * pk);

    public static BigInteger SLnqOrder(int n, int q) => GLnqOrder(n, q) / (q - 1);

    public static BigInteger GUnqOrder(int n, int q) => BigInteger.Pow(q, n * (n - 1) / 2) * n.Range(1)
        .Select(k => BigInteger.Pow(q, k) - (-1).Pow(k))
        .Aggregate(BigInteger.One, (acc, pk) => acc * pk);

    public static BigInteger SUnqOrder(int n, int q) => GUnqOrder(n, q) / (q + 1);

    public static BigInteger GOnqOrder(int n, int q)
    {
        var k = q % 2 == 0 ? 1 : 2;
        if (n % 2 != 0)
        {
            return k * BigInteger.Pow(q, (n - 1).Pow(2) / 4) *
                   ((n - 1) / 2).Range(1)
                   .Select(i => BigInteger.Pow(q, 2 * i) - 1)
                   .Aggregate(BigInteger.One, (acc, pk) => acc * pk);
        }
        else
            throw new($"n={n} must be odd");
    }

    public static BigInteger SOnqOrder(int n, int q) => q % 2 == 0 ? GOnqOrder(n, q) : GOnqOrder(n, q) / 2;

    public static BigInteger DPGLnpOrder(int n, int p) =>
        BigInteger.Pow(p - 1, n) * n.Range(1).Aggregate(BigInteger.One, (acc, k) => acc * k);

    public static BigInteger DPSLnpOrder(int n, int p) => DPGLnpOrder(n, p) / (p - 1);

    public static (MatFq a, MatFq b) GLnqGenerators(int n, int q)
    {
        var gl = new GLnq(n, q);
        var E = gl.Fq.X;
        if (n == 2)
        {
            if (q == 2)
            {
                var a = gl[1, 1, 0, 1];
                var b = gl[0, 1, 1, 0];
                return (a, b);
            }
            else
            {
                var a = gl[E, 0, 0, 1];
                var b = gl[-1, 1, -1, 0];
                return (a, b);
            }
        }
        else
        {
            var arrA = (n * n).Range().Select(_ => E.Zero).ToArray();
            var arrB = (n * n).Range().Select(_ => E.Zero).ToArray();
            arrA[0] = E;
            for (int i = 0; i < n + 1; ++i)
            {
                if (i % n != 0)
                    arrA[(n + 1) * i] = E.One;

                var k = i == 0
                    ? 0
                    : i != n
                        ? (n + 1) * i - 1
                        : n - 1;

                arrB[k] = i != n ? -E.One : E.One;
            }

            if (q == 2)
            {
                arrA[1] = E.One;
                arrB[0] = E.Zero;
            }

            var a = gl.Create(arrA);
            var b = gl.Create(arrB);
            return (a, b);
        }
    }

    public static (MatFq a, MatFq b) SLnqGenerators(int n, int q)
    {
        var gl = new GLnq(n, q);
        var E = gl.Fq.X;
        var arrA = (n * n).Range().Select(_ => E.Zero).ToArray();
        var arrB = (n * n).Range().Select(_ => E.Zero).ToArray();
        arrA[0] = E;
        arrA[n + 1] = E.Inv();
        for (int i = 0; i < n + 1; ++i)
        {
            if (i > 1 && i < n)
                arrA[(n + 1) * i] = E.One;

            var k = i == 0
                ? 0
                : i != n
                    ? (n + 1) * i - 1
                    : n - 1;

            arrB[k] = i != n ? -E.One : E.One;
        }

        if (q == 2 || q == 3)
        {
            arrA[0] = arrA[1] = arrA[n + 1] = E.One;
            arrB[0] = E.Zero;
        }

        var a = gl.Create(arrA);
        var b = gl.Create(arrB);
        return (a, b);
    }

    public static (Mat a, Mat b) GLnpGenerators(int n, int p)
    {
        var (a0, b0) = GLnqGenerators(n, p);
        var gl = new GL(n, p);
        var arrA = a0.Table.Select(e => e[0].K).ToArray();
        var arrB = b0.Table.Select(e => e[0].K).ToArray();
        return (gl.Create(arrA), gl.Create(arrB));
    }

    public static (Mat a, Mat b) SLnpGenerators(int n, int p)
    {
        var (a0, b0) = SLnqGenerators(n, p);
        var gl = new GL(n, p);
        var arrA = a0.Table.Select(e => e[0].K).ToArray();
        var arrB = b0.Table.Select(e => e[0].K).ToArray();
        return (gl.Create(arrA), gl.Create(arrB));
    }

    public const int MatrixGroupMaxOrder = 200000;

    public static ConcreteGroup<Mat> GLnp(int n, int p)
    {
        var og = GLnqOrder(n, p);
        if (og > MatrixGroupMaxOrder)
            throw new();

        var (a, b) = GLnpGenerators(n, p);
        var gl = a.GL;
        return Group.Generate(gl, a, b);
    }

    public static ConcreteGroup<Mat> SLnp(int n, int p)
    {
        var og = SLnqOrder(n, p);
        if (og > MatrixGroupMaxOrder)
            throw new();

        var (a, b) = SLnpGenerators(n, p);
        var gl = a.GL;
        return Group.Generate($"SL({n},{p})", gl, a, b);
    }

    public static ConcreteGroup<Mat> GL2p(int p) => GLnp(2, p);
    public static ConcreteGroup<Mat> SL2p(int p) => SLnp(2, p);

    public static ConcreteGroup<Coset<Mat>> L2p(int p)
    {
        if (!IntExt.Primes10000.Contains(p) || p <= 3 || p > 41)
            throw new($"p = {p} must be prime > 3 and <= 41");

        var sl2p = SL2p(p);
        var z = Group.Zentrum(sl2p);
        return sl2p.Over(z, $"L2({p})");
    }

    public static ConcreteGroup<MatFq> GLnq(int n, int q)
    {
        var og = GLnqOrder(n, q);
        if (og > MatrixGroupMaxOrder)
            throw new();

        var (a, b) = GLnqGenerators(n, q);
        var gl = a.GLnq;
        return Group.Generate(gl, a, b);
    }

    public static ConcreteGroup<MatFq> SLnq(int n, int q)
    {
        var og = SLnqOrder(n, q);
        if (og > MatrixGroupMaxOrder)
            throw new();

        var (a, b) = SLnqGenerators(n, q);
        var gl = a.GLnq;
        return Group.Generate($"SL({n},{q})", gl, a, b);
    }

    static HashSet<MatFq> GeneratorsGU2q(int q, bool special)
    {
        if (q < 2 || IntExt.PrimesDec(q).Count > 2 || (!special && GUnqOrder(2, q) > MatrixGroupMaxOrder) ||
            (special && SUnqOrder(2, q) > MatrixGroupMaxOrder))
            throw new(
                $"Out of bounds, q={q} not prime p^r or GU(2,q)>{MatrixGroupMaxOrder} or SU(2,q)>{MatrixGroupMaxOrder}");

        var q2 = q * q;
        var Glnq = new GLnq(2, q2);
        var a = Glnq.Fq.X;
        var Fq2 = Group.MulGroup($"F{q2}", a);
        var ax = a.Pow(q); // in F(q^2), (a^q)^q=a

        if (q == 2)
        {
            var J = Glnq[0, 1, 1, 0];
            if (!special)
            {
                var A = Glnq[a, a, 0, a];
                return new() { A, J };
            }
            else
            {
                var A = Glnq[1, 1, 0, 1];
                return new() { A, J };
            }
        }

        if (!special)
        {
            var gen0 = Glnq[a, 0, 0, ax.Inv()];
            var e = Fq2.Where(x => x.Equals(-x.Substitute(ax))).MinBy(x => Fq2.ElementsOrders[x]);
            var gen1 = Glnq[0, 1, 1, e];
            return new() { gen0, gen1 };
        }
        else
        {
            var e0 = Fq2.Where(x =>
            {
                var xi = x.Inv();
                var xib = xi.Substitute(ax);
                return (xi + xib).Equals(a.Zero) && (-x * xib).Equals(a.One);
            }).Distinct().MaxBy(x => Fq2.ElementsOrders[x]);
            var gen0 = Glnq[1, e0, -e0.Inv(), 0];

            if (q == 3)
            {
                var gen01 = Glnq[0, e0, -e0.Inv(), 0];
                return new() { gen0, gen01 };
            }

            var e1 = Fq2.Where(x => (x.Inv() * x.Substitute(ax)).Equals(a.One)).MaxBy(x => Fq2.ElementsOrders[x]);
            var gen1 = Glnq[e1, 0, 0, e1.Inv()];
            return new() { gen1, gen0 };
        }
    }

    static HashSet<MatFq> GeneratorsGO3q(int q, bool special)
    {
        if (q < 2 || IntExt.PrimesDec(q).Count != 1 || (!special && GOnqOrder(3, q) > MatrixGroupMaxOrder) ||
            (special && SOnqOrder(3, q) > MatrixGroupMaxOrder))
            throw new(
                $"Out of bounds, q={q} not prime p^r or GO(3,q)>{MatrixGroupMaxOrder} or SO(3,q)>{MatrixGroupMaxOrder}");

        int OrderMatOrth(EPoly<ZnInt> x0, EPoly<ZnInt> y0)
        {
            var (x1, y1) = (x0.One, y0.Zero);
            for (int i = 1; i < 50; i++)
            {
                (x1, y1) = (x0 * x1 - y0 * y1, x0 * y1 + y0 * x1);
                if (x1.Equals(x1.One) && y1.IsZero())
                    return i;
            }

            throw new("####################################");
        }

        var Glnq = new GLnq(3, q);
        var a = Glnq.Fq.X;
        var arrFq = Group.MulGroup($"F{q}", a);
        var pows = q.Range(-1).ToDictionary(k => k == -1 ? a.Zero : a.Pow(k), k => k);
        var square = arrFq.Append(a.Zero).Select(x => (x, x2: x * x)).GroupBy(e => e.x2)
            .ToDictionary(e => e.Key, e => e.Select(f => f.x).ToArray());
        var dicSquare =
            arrFq.ToDictionary(x => x, x => square.TryGetValue(x, out var value) ? value : new EPoly<ZnInt>[0]);
        dicSquare[a.Zero] = new EPoly<ZnInt>[0];
        var possibles = arrFq.Append(a.Zero)
            .Select(x0 => (x: x0, yList: dicSquare[1 - x0 * x0]))
            .Where(e => e.yList.Length != 0)
            .SelectMany(e => e.yList.Select(y0 => (e.x, y: y0)))
            .Distinct()
            .Select(e => (e.x, e.y, mat: Glnq[1, 0, 0, 0, e.x, e.y, 0, -e.y, e.x]))
            .Select(e => (e, OrderMatOrth(e.x, e.y)))
            .OrderByDescending(e => e.Item2)
            .ThenByDescending(e => pows[e.e.x])
            .ToArray();

        var e0 = special ? a.One : -a.One;
        var ide = Glnq[e0, 0, 0, 0, e0, 0, 0, 0, e0];

        var m0 = possibles[0].e.mat;
        var m1 = q != 5 ? Glnq[0, 1, 0, 0, 0, 1, 1, 0, 0] : Glnq[3, 1, 1, 1, 4, 3, 4, 2, 1];
        return new() { m0, Glnq.Op(m1, ide) };
    }

    public static ConcreteGroup<MatFq> GU2q(int q)
    {
        var gens = GeneratorsGU2q(q, special: false);
        var Glnq = gens.First().GLnq;
        return Group.Generate($"GU(2,{q})", Glnq, [..gens]);
    }

    public static ConcreteGroup<MatFq> GO3q(int q)
    {
        var gens = GeneratorsGO3q(q, special: false);
        var Glnq = gens.First().GLnq;
        return Group.Generate($"GO(3,{q})", Glnq, [..gens]);
    }

    public static ConcreteGroup<Mat> GO3p(int p)
    {
        if (!IntExt.Primes10000.Contains(p))
            throw new($"p={p} isnt prime");

        var gens = GeneratorsGO3q(p, special: false);
        var gl = new GL(3, p);
        var gens0 = gens.Select(m => gl.Create(m.Table.Select(c => c[0].K).ToArray())).ToArray();
        return Group.Generate($"GO(3,{p})", gl, gens0);
    }

    public static ConcreteGroup<Mat> SO3p(int p)
    {
        if (!IntExt.Primes10000.Contains(p))
            throw new($"p={p} isnt prime");

        var gens = GeneratorsGO3q(p, special: true);
        var gl = new GL(3, p);
        var gens0 = gens.Select(m => gl.Create(m.Table.Select(c => c[0].K).ToArray())).ToArray();
        return Group.Generate($"SO(3,{p})", gl, gens0);
    }

    public static ConcreteGroup<MatFq> SU2q(int q)
    {
        var gens = GeneratorsGU2q(q, special: true);
        var Glnq = gens.First().GLnq;
        return Group.Generate($"SU(2,{q})", Glnq, [..gens]);
    }

    public static ConcreteGroup<MatFq> SO3q(int q)
    {
        var gens = GeneratorsGO3q(q, special: true);
        var Glnq = gens.First().GLnq;
        return Group.Generate($"SO(3,{q})", Glnq, [..gens]);
    }

    public static Mat[] DPGLnpGenerators(int n, int p)
    {
        var og = DPGLnpOrder(n, p);
        if (og > MatrixGroupMaxOrder)
            throw new();

        if (!IntExt.Primes10000.Contains(p))
            throw new($"p={p} isnt prime");


        var gl = new GL(n, p);
        var e0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(gl.P, p - 1);
        var diag = gl.At(gl.Neutral().Table, 0, e0);
        return SnGensMat(n).Select(e => gl.Create(e.Table)).Prepend(diag).ToArray();
    }

    public static Mat[] DPSLnpGenerators(int n, int p)
    {
        var og = DPSLnpOrder(n, p);
        if (og > MatrixGroupMaxOrder)
            throw new();

        if (!IntExt.Primes10000.Contains(p))
            throw new($"p={p} isnt prime");

        var gl = new GL(n, p);
        var e0 = IntExt.Solve_k_pow_m_equal_one_mod_n_strict(gl.P, p - 1);
        var e1 = IntExt.InvModPbez(e0, p);
        var diag = gl.At(gl.Neutral().Table, (0, n + 1), (e0, e1));
        return SnGensMat(n).Select(e => gl.Create(e.Table)).Select(m =>
        {
            var det = m.Det;
            var m0 = gl.At(gl.Neutral().Table, 0, IntExt.InvModPbez(det, p));
            return gl.Op(m0, m);
        }).Prepend(diag).ToArray();
    }

    public static MatFq[] DPGLnqGenerators(int n, int q)
    {
        var og = FG.DPGLnpOrder(n, q);
        if (og > FG.MatrixGroupMaxOrder || IntExt.PrimesDecomposition(q).Distinct().Count() != 1)
            throw new();

        var gl = new GLnq(n, q);
        Console.WriteLine(gl);
        var x = gl.Fq.X;
        var fq = Group.MulGroup("Fq", x);
        var e0 = fq.First(e => fq.ElementsOrders[e] == q - 1);
        var diag = gl.Neutral().Table.Select((c, i) => i == 0 ? e0 : c).ToArray();
        return FG.SnGensMat(n).Select(e => gl.Create(e.Table.Select(k => k * x.One).ToArray()))
            .Prepend(gl.Create(diag)).ToArray();
    }

    public static MatFq[] DPSLnqGenerators(int n, int q)
    {
        var og = FG.DPSLnpOrder(n, q);
        if (og > FG.MatrixGroupMaxOrder)
            throw new();

        var gl = new GLnq(n, q);
        var x = gl.Fq.X;
        var fq = Group.MulGroup("Fq", x);
        var e0 = fq.First(e => fq.ElementsOrders[e] == q - 1);
        var e1 = e0.Inv();
        var diag = gl.Neutral().Table.Select((c, i) => i == 0 ? e0 : (i == n + 1 ? e1 : c)).ToArray();
        return FG.SnGensMat(n).Select(e => gl.Create(e.Table.Select(k => k * x.One).ToArray())).Select(m =>
        {
            var deti = m.Det.Inv();
            var m0 = gl.Create(gl.Neutral().Table.Select((c, i) => i == 0 ? deti : c).ToArray());
            return gl.Op(m0, m);
        }).Prepend(gl.Create(diag)).ToArray();
    }

    public static ConcreteGroup<Mat> DPGLnp(int n, int p)
    {
        var gens = DPGLnpGenerators(n, p);
        var gl = gens[0].GL;
        return Group.Generate($"DPGL({n},{p})", gl, gens);
    }

    public static ConcreteGroup<Mat> DPSLnp(int n, int p)
    {
        var gens = DPSLnpGenerators(n, p);
        var gl = gens[0].GL;
        return Group.Generate($"DPSL({n},{p})", gl, gens);
    }

    public static ConcreteGroup<MatFq> DPGLnq(int n, int q)
    {
        var gens = DPGLnqGenerators(n, q);
        var gl = gens[0].GLnq;
        return Group.Generate($"DPGL({n},F{q})", gl, gens);
    }

    public static ConcreteGroup<MatFq> DPSLnq(int n, int q)
    {
        var gens = DPSLnqGenerators(n, q);
        var gl = gens[0].GLnq;
        return Group.Generate($"DPSL({n},F{q})", gl, gens);
    }

    public static ConcreteGroup<Mat> AbelianMat(params int[] seq)
    {
        var dim = seq.Length;
        var p = IntExt.Primes10000.First(p => seq.All(o => (p - 1) % o == 0));
        var gl = new GL(dim, p);
        var Up = FG.UnInt(p);
        var seq2 = seq.Select(o => Up.ElementsOrders.First(e => e.Value == o).Key.K).ToArray();
        var gens = seq2.Select((v, k) => gl.At(gl.Neutral().Table, k * (dim + 1), v)).ToArray();
        return Group.Generate(seq.Glue(" x ", "C{0}"), gl, gens);
    }

    public static Mat[] SnGensMat(int n, int p = 2)
    {
        if (n < 1 || !IntExt.Primes10000.Contains(p))
            throw new GroupException(GroupExceptionType.GroupDef);

        var sn = new Sn(n);
        var gl = new GL(n, p);
        var id = gl.Neutral().Table;
        var idc = id.Chunk(gl.N).ToArray();
        return sn.GetGenerators().Select(e => gl.Create(e.Apply(idc).SelectMany(l => l).ToArray())).ToArray();
    }

    public static Mat[] AnGensMat(int n, int p = 2)
    {
        if (n < 3 || !IntExt.Primes10000.Contains(p))
            throw new GroupException(GroupExceptionType.GroupDef);

        var sn = new Sn(n);
        var gl = new GL(n, p);
        var id = gl.Neutral().Table;
        var idc = id.Chunk(gl.N).ToArray();
        var gensAn = (gl.N - 2).Range(3).Select(i => sn[(1, 2, i)]).ToArray();
        return gensAn.Select(e => gl.Create(e.Apply(idc).SelectMany(l => l).ToArray())).ToArray();
    }

    public static ConcreteGroup<Mat> SymmetricMat(int n)
    {
        var gensMatSn = SnGensMat(n);
        var gl = gensMatSn[0].GL;
        return Group.Generate($"Symm{n}", gl, gensMatSn);
    }

    public static ConcreteGroup<Mat> AlternateMat(int n)
    {
        var gensMatAn = AnGensMat(n);
        var gl = gensMatAn[0].GL;
        return Group.Generate($"Alt{n}", gl, gensMatAn);
    }

    public static ConcreteGroup<Mat> MetaCyclicGLnp_DiagByPerm(int m, int n, int r, int dim)
    {
        // Console.WriteLine($"Solve M({m}x:{n}){r} in GL({dim},{p})");
        var distinctTypes = IntExt.Partitions32[dim].Select(l => l.Order().ToArray()).OrderBy(l => l.Length).ToArray();
        var nks = distinctTypes.Select(l => l.Aggregate((a0, a1) => a0 * a1))
            .SelectMany(e => IntExt.Dividors(e).Append(e).Where(j => j != 1)).Append(n).ToHashSet();
        foreach (var p in nks.Select(nk => IntExt.Primes10000.First(p => (p - 1) % m == 0 && (p - 1) % nk == 0))
                     .Distinct().Order())
        {
            var Up = UnInt(p);
            var gl = new GL(dim, p);
            var ordn = Up.Where(e => n % Up.ElementsOrders[e] == 0)
                .OrderBy(e => Up.ElementsOrders[e])
                .Select(e => gl.At(gl.Neutral().Table, 0, e.K))
                .ToArray();
            var sn = new Sn(dim);
            var m1s = IntExt.Partitions32[dim].OrderBy(l => l.Count)
                .Select(t => IntExt.PermAndCyclesFromType(t.Order().ToArray()))
                .Select<(int[] perm, int[][] cycles), IEnumerable<(int[] perm, int[][] cycles, Mat mat)>>(e =>
                {
                    var e0 = gl.Neutral().Table.Chunk(dim).ToArray();
                    var perm = sn.CreateElement(e.perm.Select(i => i + 1).ToArray());
                    var e1 = perm.Apply(e0);
                    var mat0 = gl.Create(e1.SelectMany(v => v).ToArray());
                    return ordn.Select(mat => gl.Op(mat0, mat))
                        .Where(mat => mat.IsOrder(n))
                        .Select(mat => (e.perm, e.cycles, mat));
                })
                .SelectMany(e => e);

            foreach (var (perm, cycles, m1) in m1s)
            {
                var m1i = gl.Invert(m1);
                var seq = cycles.Select(c => c.Length).Select(l =>
                {
                    var r0 = IntExt.PowMod(r, l, m);
                    var ordm = Up.Where(e => m % Up.ElementsOrders[e] == 0 && e.Pow(r0).Equals(e))
                        .OrderByDescending(e => Up.ElementsOrders[e]);
                    return ordm.Select(a => l.Range(1).Select(k => a.Pow(IntExt.PowMod(r, k, m))).Reverse().ToArray());
                }).MultiLoop().Select(l => l.ToArray());
                foreach (var l in seq)
                {
                    var arr = new int[dim * dim];
                    foreach (var (sols, idxs) in l.Zip(cycles))
                    foreach (var (idx, sol) in idxs.Zip(sols))
                        arr[perm[idx] * (dim + 1)] = sol.K;

                    var m0 = gl.Create(arr);
                    if (m0.IsOrder(m) && gl.Op(m1i, gl.Op(m0, m1)).Equals(gl.Times(m0, r)))
                    {
                        var name = IntExt.Gcd(m, n * (r - 1)) == 1 ? $"F({m}x:{n}){r}" : $"M({m}x:{n}){r}";
                        var mtGL = Group.Generate(name, gl, m0, m1);
                        if (mtGL.Count() == m * n)
                            return mtGL;
                    }
                }
            }
        }

        return Group.Generate(new GL(1, 2));
    }

    public static ConcreteGroup<Mat> MetaCyclicSdpMat(int m, int n, int r, int maxDim = 12)
    {
        foreach (var dim in maxDim.Range(1).Where(d => d != 1 && (IntExt.Gcd(m, d) != 1 || IntExt.Gcd(m - 1, d) != 1)))
        {
            var mtGL = MetaCyclicGLnp_DiagByPerm(m, n, r, dim);
            if (mtGL.Count() != 1)
                return mtGL;
        }

        throw new GroupException(GroupExceptionType.GroupDef);
    }

    public static List<ConcreteGroup<Mat>> MetaCyclicSdpMat(int ord, int maxDim = 12)
    {
        return IntExt.Dividors(ord).Where(d => d > 1)
            .SelectMany(m => MetaCyclicSdpGetR(m, ord / m).Select(r => (m, n: ord / m, r)))
            .Select(e => MetaCyclicSdpMat(e.m, e.n, e.r, maxDim))
            .ToList();
    }

    public static GLn<K> GLnK<K>(int n, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new GLn<K>(n, scalar);
    }

    public static GLn<K> GLnK<K>(string name, int n, K scalar) where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    {
        return new GLn<K>(name, n, scalar);
    }

    public static Polynomial<ZnInt, T> ToZnInt<T>(this Polynomial<ZnInt, T> P, int mod) where T : struct, IElt<T>
    {
        return new Polynomial<ZnInt, T>(P.Indeterminates, ZnInt.ZnZero(mod),
            new(P.Coefs.ToDictionary(kv => kv.Key, kv => new ZnInt(mod, kv.Value.K))));
    }

    public static Polynomial<ZnBInt, T> ToZnBInt<T>(this Polynomial<ZnInt, T> P) where T : struct, IElt<T>
    {
        var mod = P.KZero.Mod;
        return new Polynomial<ZnBInt, T>(P.Indeterminates, ZnBInt.ZnZero(mod),
            new(P.Coefs.ToDictionary(kv => kv.Key, kv => new ZnBInt(mod, kv.Value.K))));
    }

    public static Polynomial<ZnInt, T> ToZnInt<T>(this Polynomial<ZnBInt, T> P) where T : struct, IElt<T>
    {
        if (P.KZero.Details.O > 1)
            throw new();
        var mod = (int)P.KZero.Mod;
        return new Polynomial<ZnInt, T>(P.Indeterminates, ZnInt.ZnZero(mod),
            new(P.Coefs.ToDictionary(kv => kv.Key, kv => new ZnInt(mod, (int)kv.Value.K))));
    }

    public static Polynomial<Rational, T> ToRationalPoly<T>(this Polynomial<ZnInt, T> P) where T : struct, IElt<T>
    {
        return new Polynomial<Rational, T>(P.Indeterminates, Rational.KZero(),
            new(P.Coefs.ToDictionary(kv => kv.Key, kv => new Rational(kv.Value.K))));
    }

    public static KMatrix<ZnInt> ToKMatrix(this Mat m) => m.Table.Select(e => new ZnInt(m.GL.P, e)).ToKMatrix(m.GL.N);

    public static GLn<ZnInt> ToGLnK(this GL gl) => GLnK($"F{gl.P}", gl.N, ZnInt.ZnZero(gl.P));

    public static Polynomial<ZnInt, Xi> Mod(this Polynomial<ZnInt, Xi> P, int mod)
    {
        var coefs = P.Coefs.Select(e => (e.Key, IntExt.AmodP(e.Value.K, mod)))
            .Where(e => e.Item2 != 0)
            .ToDictionary(e => e.Key, e => new ZnInt(P.P, e.Item2));
        return new(P.Indeterminates, new ZnInt(P.P, 0), new(coefs));
    }

    public static EPolynomial<Rational>[] NumberFieldQ(Polynomial<Rational, Xi>[] basis)
    {
        var e0 = basis[0];
        var polyBasis = new PolynomialBasis<Rational, Xi>(e0.Indeterminates, basis);
        return e0.Indeterminates
            .Select(xi => new Polynomial<Rational, Xi>(new Monom<Xi>(e0.Indeterminates, xi, 1), e0.KOne))
            .Select(xi => new EPolynomial<Rational>(xi, polyBasis)).ToArray();
    }

    public static EPolynomial<Rational> NumberFieldQ(Polynomial<Rational, Xi> e0) => NumberFieldQ(new[] { e0 })[0];

    public static (EPolynomial<Rational>, EPolynomial<Rational>) NumberFieldQ(Polynomial<Rational, Xi> e0,
        Polynomial<Rational, Xi> e1)
    {
        var nbf = NumberFieldQ(new[] { e0, e1 });
        var x0 = e0.ExtractIndeterminate;
        var i0 = e0.Indeterminates.ToList().FindIndex(xi => xi.Equals(x0));
        var x1 = e1.ExtractIndeterminate;
        var i1 = e0.Indeterminates.ToList().FindIndex(xi => xi.Equals(x1));
        return (nbf[i0], nbf[i1]);
    }

    public static EPolynomial<Rational> NumberFieldQ(KPoly<Rational> e, string x)
    {
        var (_, t) = Ring.Polynomial(Rational.KZero(), x, "_t_").Deconstruct();
        var a = e.ToPolynomial(t.Indeterminates, t.Indeterminates[0]);
        return NumberFieldQ(a);
    }

    public static (EPolynomial<Rational> x, EPolynomial<Rational> p0) NumberFieldQ(KPoly<Rational> e, string x,
        string p0)
    {
        var all = new[] { x, p0, "_t_" };
        var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, all);
        var x0 = xis[0];
        var a = e.ToPolynomial(x0.Indeterminates, x0.Indeterminates[0]);
        var nf = NumberFieldQ(a);
        return (nf, new EPolynomial<Rational>(xis[1], nf.Basis));
    }

    public static (EPolynomial<Rational>, EPolynomial<Rational>) NumberFieldQ((KPoly<Rational>, string) e0,
        (KPoly<Rational>, string) e1)
    {
        var (_, _, t) = Ring.Polynomial(Rational.KZero(), e0.Item2, e1.Item2, "_t_").Deconstruct();
        var a = e0.Item1.ToPolynomial(t.Indeterminates, t.Indeterminates[0]);
        var b = e1.Item1.ToPolynomial(t.Indeterminates, t.Indeterminates[1]);
        var nbf = NumberFieldQ(new[] { a, b });
        return (nbf[0], nbf[1]);
    }

    public static (EPolynomial<Rational> e0, EPolynomial<Rational> e1, EPolynomial<Rational> p0) NumberFieldQ(
        (KPoly<Rational>, string) e0,
        (KPoly<Rational>, string) e1, string p0, params string[] others)
    {
        var all = new[] { e0.Item2, e1.Item2, p0, "_t_" }.ToArray();
        var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, all);
        var x0 = xis[0];
        var a = e0.Item1.ToPolynomial(x0.Indeterminates, x0.Indeterminates[0]);
        var b = e1.Item1.ToPolynomial(x0.Indeterminates, x0.Indeterminates[1]);
        var nbf = NumberFieldQ(new[] { a, b });
        return (nbf[0], nbf[1], new EPolynomial<Rational>(xis[2], nbf[0].Basis));
    }

    // Algebre Tome 1, page 415
    public static KPoly<Rational>[] GaloisGroupPolynomialsList(int n)
    {
        var x = QPoly();
        if (n == 1)
            return new[] { x };
        if (n == 2)
            return new[] { x.Pow(2) + x + 1 };
        if (n == 3)
            return new[]
            {
                x.Pow(3) + x.Pow(2) - 2 * x - 1,
                x.Pow(3) + 2
            };
        if (n == 4)
            return new[]
            {
                x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1,
                x.Pow(4) + 1, x.Pow(4) - 2,
                x.Pow(4) + 8 * x + 12,
                x.Pow(4) + x + 1
            };
        if (n == 5)
            return new[]
            {
                x.Pow(5) + x.Pow(4) - 4 * x.Pow(3) - 3 * x.Pow(2) + 3 * x + 1,
                x.Pow(5) - 5 * x + 12,
                x.Pow(5) + 2,
                x.Pow(5) + 20 * x + 16,
                x.Pow(5) - x + 1
            };
        if (n == 6)
            return new[]
            {
                x.Pow(6) + x.Pow(5) + x.Pow(4) + x.Pow(3) + x.Pow(2) + x + 1, x.Pow(6) + 108, x.Pow(6) + 2,
                x.Pow(6) - 3 * x.Pow(2) - 1, x.Pow(6) + 3 * x.Pow(3) + 3,
                x.Pow(6) - 3 * x.Pow(2) + 1, x.Pow(6) - 4 * x.Pow(2) - 1,
                x.Pow(6) - 3 * x.Pow(5) + 6 * x.Pow(4) - 7 * x.Pow(3) + 2 * x.Pow(2) + x - 4,
                x.Pow(6) + 2 * x.Pow(3) - 2,
                x.Pow(6) + 6 * x.Pow(4) + 2 * x.Pow(3) + 9 * x.Pow(2) + 6 * x - 4,
                x.Pow(6) + 2 * x.Pow(2) + 2,
                x.Pow(6) - 2 * x.Pow(5) - 5 * x.Pow(2) - 2 * x - 1,
                x.Pow(6) + 2 * x.Pow(4) + 2 * x.Pow(3) + x.Pow(2) + 2 * x + 2,
                x.Pow(6) - x.Pow(5) - 10 * x.Pow(4) + 30 * x.Pow(3) - 31 * x.Pow(2) + 7 * x + 9,
                x.Pow(6) + 24 * x - 20,
                x.Pow(6) + x + 1
            };
        if (n == 7)
            return new[]
            {
                x.Pow(7) + x.Pow(6) - 12 * x.Pow(5) - 7 * x.Pow(4) + 28 * x.Pow(3) + 14 * x.Pow(2) - 9 * x + 1,
                x.Pow(7) + 7 * x.Pow(3) + 7 * x.Pow(2) + 7 * x - 1,
                x.Pow(7) - 14 * x.Pow(5) + 56 * x.Pow(3) - 56 * x + 22,
                x.Pow(7) + 2,
                x.Pow(7) - 7 * x.Pow(3) + 14 * x.Pow(2) - 7 * x + 1, x.Pow(7) + 7 * x.Pow(4) + 14 * x + 3,
                x.Pow(7) + x + 1
            };

        throw new();
    }

}