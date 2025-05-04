using System.Numerics;
using System.Text;
using System.Text.RegularExpressions;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;
using Group = FastGoat.Structures.Group;
using RegXGroup = System.Text.RegularExpressions.Group;
using GFelt = FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>;
using Rl = FastGoat.Structures.VecSpace.EPoly<FastGoat.Structures.VecSpace.EPoly<FastGoat.UserGroup.Integers.ZnInt>>;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RecomputeAllPrimesUpTo(200000);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;
var setGF = new HashSet<(int p, int n)>();
var (nbErrors, nbAp, nbCurves) = (0, 0, 0);
var Infos = new StringBuilder();

Dictionary<int, List<int>> SmallPrimesList(Rational N)
{
    var nmax = 2 * double.Sqrt(N);
    var listL = new Dictionary<int, List<int>>();
    foreach (var p in Primes10000.Where(p => p <= nmax && !N.Mod(p).IsZero()))
    {
        var lmax = 4 * double.Sqrt(p);
        listL[p] = new();
        var L = 1;
        foreach (var l in Primes10000.Where(l => p % l != 0))
        {
            listL[p].Add(l);
            L *= l;
            if (L > lmax)
                break;
        }
    }

    return listL;
}

EllPt<GFelt> Frob(EllPt<GFelt> pt, int n = 1)
{
    if (pt.IsO)
        return pt;

    var p = pt.X.P;
    return new(pt.X.FastPow(BigInteger.Pow(p, n)), pt.Y.FastPow(BigInteger.Pow(p, n)));
}

EllPt<GFelt> FrobEq(EllGroup<GFelt> E, EllPt<GFelt> pt, int t)
{
    if (pt.IsO)
        return pt;

    var p = pt.X.P;
    var phi2 = Frob(Frob(pt));
    var t_phi = Frob(E.Times(pt, -t));
    var q = E.Times(pt, p);
    return E.Op(phi2, E.Op(t_phi, q));
}

bool Relation(EllGroup<GFelt> E, EllPt<GFelt> pt, int t)
{
    if (pt.IsO)
        return true;

    var p = pt.X.P;
    var phi2 = Frob(Frob(pt));
    var t_phi = Frob(E.Times(pt, -t));
    var q = E.Times(pt, p);
    return E.Op(phi2, E.Op(t_phi, q)).IsO;
}

FracPoly<FracPoly<ZnInt>> FastPowRl(FracPoly<FracPoly<ZnInt>> a, BigInteger k, KPoly<ZnInt> R0, KPoly<ZnInt> R1,
    KPoly<FracPoly<ZnInt>> R2, KPoly<FracPoly<ZnInt>> R3)
{
    if (k == 0)
        return a.One;

    if (k < 0)
        return FastPowRl(a.Inv(), -k, R0, R1, R2, R3);

    var (r, a0, e0) = (a.One, a, k);
    while (e0 > 0)
    {
        if (e0 % 2 == 1)
            r = Simplify5(r * a0, R0, R1, R2, R3);

        e0 >>= 1;
        a0 = Simplify5(a0 * a0, R0, R1, R2, R3);
    }

    return r;
}

EllPt<FracPoly<FracPoly<ZnInt>>> FrobRl(EllPt<FracPoly<FracPoly<ZnInt>>> pt, int n, KPoly<ZnInt> R0, KPoly<ZnInt> R1,
    KPoly<FracPoly<ZnInt>> R2, KPoly<FracPoly<ZnInt>> R3)
{
    if (pt.IsO)
        return pt;

    var p = pt.X.P;
    return new(FastPowRl(pt.X, BigInteger.Pow(p, n), R0, R1, R2, R3),
        FastPowRl(pt.Y, BigInteger.Pow(p, n), R0, R1, R2, R3));
}

int FrobTrace(EllGroup<GFelt> E, EllPt<GFelt>[] nTors, int l) =>
    l.SeqLazy().First(t => nTors.All(pt => Relation(E, pt, t)));

void EllApFrobTrace(BigInteger[] curve)
{
    GlobalStopWatch.AddLap();
    var E = EC.EllCoefs(curve);
    var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
    var Ell = E.ToEllGroup();
    var allList = SmallPrimesList(N);
    var pmax = allList.Keys.Max();
    var lmax = allList.Max(e => e.Value.Max());
    Console.WriteLine($"{Ell} {Ell.ShortFormStr}");

    var divPolys = EC.DivisionPolynomial(Ell, lmax + 3).f.ToDictionary(e => e.Key, e => e.Value.PrimitiveZPoly());
    divPolys.Println("divPolys");

    Console.WriteLine($"N = {N} pmax = {pmax} listMax = {lmax}");
    foreach (var (p, listL) in allList)
    {
        if (p < 5)
            continue;

        var ap = EC.EllAp(Ell, p);
        var frobTr = new Dictionary<int, int>();
        foreach (var l in listL)
        {
            if (p % l == 0)
                continue;

            var Ep = Ell.ToGF(p);
            Ep.Field = $"GF({p})";
            Console.WriteLine($"{l}-torsion of {Ep} from {Ell}");
            Console.WriteLine($"{l}-divPol mod {p} = {divPolys[l].ToZnPoly(p)}");
            var (nTors, g) = EC.EllFpNTors(Ell, divPolys[l], p, l);
            l.SeqLazy().Grid2D(nTors).Select(e => (e, FrobEq(Ep.ToGF(g.X), e.t2, e.t1))).Println();
            var Egf = Ell.ToGF(g);
            setGF.Add((p, g.X.F.Degree));
            frobTr[l] = FrobTrace(Egf, nTors.GetGenerators().ToArray(), l);
        }

        frobTr.Println("Frob Traces");
        var keys = frobTr.Keys.ToArray();
        var values = keys.Select(k => frobTr[k]).ToArray();
        var crtTable = NumberTheory.CrtTable(keys);
        var L = keys.Aggregate((li, lj) => li * lj);
        var crt = NumberTheory.CRT(values, crtTable, L);
        var ap1 = crt < L / 2 ? crt : crt - L;
        Console.WriteLine($"p = {p} ap = {ap} crt = {crt} ap1 = {ap1} L = {L} Check:{ap == ap1}");
        if (ap != ap1)
        {
            ++nbErrors;
            Console.WriteLine($"Warning[{nbErrors}]");
        }

        ++nbAp;
        Console.WriteLine();
    }

    ++nbCurves;
    GlobalStopWatch.Show($"End FrobTrace {Ell}");
    Console.WriteLine();
}

void testEllApFrobTrace1()
{
    EllApFrobTrace([1, 0]);
    EllApFrobTrace([-1, 0]);
    EllApFrobTrace([1, 1]);
    EllApFrobTrace([-43, 166]);
}

void testEllApFrobTrace2()
{
    EllApFrobTrace([1, 0]);
    EllApFrobTrace([-1, 0]);
    EllApFrobTrace([1, 1]);
    EllApFrobTrace([-43, 166]);
    EllApFrobTrace([0, 0, 1, 0, -7]);
    EllApFrobTrace([1, 0, 0, 0, 1]);
    EllApFrobTrace([1, -1, 0, -4, 4]);
    EllApFrobTrace([1, -1, 1, -19353, 958713]);
    EllApFrobTrace([1, 1, 1, -17714, 900047]);
}

void testEllApFrobTrace3()
{
    foreach (var e in EllipticExt.LMFDB_Ell_Q().Where(e => e.conductor < 1000))
        EllApFrobTrace(e.model);
}

void runTest(Action act)
{
    GlobalStopWatch.AddLap();
    (nbErrors, nbAp, nbCurves) = (0, 0, 0);
    setGF.Clear();
    act();
    var (p, n) = setGF.MaxBy(e => BigInteger.Pow(e.p, e.n));
    Console.WriteLine($"{act.Method.Name}");
    Console.WriteLine($"Nb Curves {nbCurves}");
    Console.WriteLine($"Max field GF({p}^{n})");
    Console.WriteLine($"Total Errors {nbErrors}/{nbAp}");
    GlobalStopWatch.Show();
    setGF.Where(e => !PolynomExt.AllConwayPolys.ContainsKey(e.p) || !PolynomExt.AllConwayPolys[e.p].ContainsKey(e.n))
        .Println(gf => $"GF({gf.p}^{gf.n})", "Missing Fq");
    Console.WriteLine();
}

void runAll()
{
    GlobalStopWatch.Restart();
    runTest(testEllApFrobTrace1);
    runTest(testEllApFrobTrace2);
    runTest(testEllApFrobTrace3);

    // 5-torsion of Ell[1,1](GF(13)) from Ell[1,1](Q)
    // 5-divPol mod 13 =  5*X^12 + 10*X^10 +  3*X^9 + 12*X^8 +  6*X^7 +  6*X^6 +  6*X^5 +  9*X^4 + 10*X^3 +  9*X^2 + X + 12
    // x(P) and y(P) in GF(13^4)
    //     MinPol   of GF(13^4) = X^4 + 11*X^3 +  4*X^2 +  6*X +  7
    //     PrimRoot of GF(13^4) =  7*a^3 + 11*a^2 +  6*a
    // |5-Tors(Ell[1,1](GF(13^4)))| = 25
    // Type        AbelianGroup
    // BaseGroup   Ell[1,1](GF(13^4))
    // 
    // Generators of 5-Tors(Ell[1,1](GF(13^4)))
    // gen1 of order 5
    // ( 2, 11*a^3 + 10*a^2 + 10*a + 1)
    // gen2 of order 5
    // ( 6*a^3 +  9*a^2 +  9*a +  8,  6*a^3 + 11*a^2 +  8*a +  6)
    // 
    // 5-Tors(Ell[1,1](GF(13^4))) ~ C5 x C5
    //
    // ...
    //
    // testEllApFrobTrace1
    // Nb Curves 4
    // Max field GF(5^48)
    // Total Errors 0/20
    // #  Time:1m28s
    // 
    // testEllApFrobTrace2
    // Nb Curves 9
    // Max field GF(61^48)
    // Total Errors 0/69
    // #  Time:4m24s
    // 
    // testEllApFrobTrace3
    // Nb Curves 88
    // Max field GF(61^48)
    // Total Errors 0/427
    // #  Time:19m30s
    // 
}

KPoly<ZnInt> Simplify1(KPoly<ZnInt> P, KPoly<ZnInt> R0, KPoly<ZnInt> R1)
{
    var (quo, rem) = P.Div(R0);
    if (quo.IsZero())
        return P;

    return Simplify1(quo * R1 + rem, R0, R1);
}

KPoly<FracPoly<ZnInt>> Simplify2(KPoly<FracPoly<ZnInt>> P, KPoly<FracPoly<ZnInt>> R0, KPoly<FracPoly<ZnInt>> R1)
{
    var (quo, rem) = P.Div(R0);
    if (quo.IsZero())
        return P;

    return Simplify2(quo * R1 + rem, R0, R1);
}

FracPoly<ZnInt> Simplify3(FracPoly<ZnInt> P, KPoly<ZnInt> R0, KPoly<ZnInt> R1)
{
    var num = Simplify1(P.Num, R0, R1);
    var denom = Simplify1(P.Denom, R0, R1);
    return new(num, denom);
}

FracPoly<FracPoly<ZnInt>> Simplify4(FracPoly<FracPoly<ZnInt>> P, KPoly<ZnInt> R0, KPoly<ZnInt> R1)
{
    var num = new KPoly<FracPoly<ZnInt>>(P.x, P.Num.KZero, P.Num.Coefs.Select(n0 => Simplify3(n0, R0, R1)).ToArray());
    var denom = new KPoly<FracPoly<ZnInt>>(P.x, P.Denom.KZero,
        P.Denom.Coefs.Select(n0 => Simplify3(n0, R0, R1)).ToArray());
    return new(num, denom);
}

FracPoly<FracPoly<ZnInt>> Simplify5(FracPoly<FracPoly<ZnInt>> P, KPoly<ZnInt> R0, KPoly<ZnInt> R1,
    KPoly<FracPoly<ZnInt>> R2, KPoly<FracPoly<ZnInt>> R3)
{
    // Console.WriteLine(new { P });
    var P1 = Simplify4(P, R0, R1);
    // Console.WriteLine(new { P1 });
    var num = Simplify2(P1.Num, R2, R3);
    var denom = Simplify2(P1.Denom, R2, R3);
    var P2 = new FracPoly<FracPoly<ZnInt>>(num, denom);
    // Console.WriteLine(new { P2 });
    var P3 = Simplify4(P2, R0, R1);
    // Console.WriteLine(new { P3 });
    return P3;
}

void RingRl(BigInteger[] curve)
{
    var E = EC.EllCoefs(curve);
    var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
    var Ell = E.ToEllGroup();
    var allList = SmallPrimesList(N);
    var pmax = allList.Keys.Max();
    var lmax = allList.Max(e => e.Value.Max());
    Console.WriteLine($"{Ell} {Ell.ShortFormStr}");

    var divPolys = EC.DivisionPolynomial(Ell, lmax + 3).f.ToDictionary(e => e.Key, e => e.Value.PrimitiveZPoly());
    divPolys.Println("divPolys");

    Console.WriteLine($"N = {N} pmax = {pmax} listMax = {lmax}");
    foreach (var (p, listL) in allList)
    {
        if (p < 5)
            continue;

        var a0 = NumberTheory.PrimitiveRootMod(p);
        var ap = EC.EllAp(Ell, p);
        var frobTr = new Dictionary<int, int>();
        foreach (var l in listL)
        {
            Console.WriteLine($"p={p} l={l}");
            var psi = divPolys[l].ToZnPoly(p).SubstituteChar('X');
            Console.WriteLine($"psi = {psi}");
            var facts = IntFactorisation.MusserSFF(psi)
                .SelectMany(e => IntFactorisation.CantorZassenhausAECF(e.g, new(p, a0), p))
                .OrderBy(f => f.Degree).ToArray();
            Console.WriteLine($"psi ~ {facts.Glue(" x ", "({0})")}");
            var f = facts.Last();
            var x = FG.KFracPoly(f.X).X;
            var Y = FG.KFracPoly('Y', x).X;
            var X = Y.KOne.X + Y.Zero;
            var Psi = psi.Substitute(Y.KOne.X) + Y.Zero;
            Console.WriteLine(Psi);
            var (a1, a2, a3, a4, a6) = E.ToEllGroup().ToZnInt(p).ArrCoefs.Select(c => c * Y.KOne * Y.One).Deconstruct();
            var lhs = Y.Pow(2) + a1 * X * Y + a3 * Y;
            var rhs = Simplify4(X.Pow(3) + a2 * X.Pow(2) + a4 * X + a6, psi, psi.Zero);
            var Erl = new EllGroup<FracPoly<FracPoly<ZnInt>>>(a1, a2, a3, a4, a6);
            Console.WriteLine(new { lhs, rhs });
            if (rhs.IsZero())
            {
                frobTr[l] = 0;
                continue;
            }
            Erl.Simplify = P0 => Simplify5(P0, psi, psi.Zero, lhs.Num, rhs.Num);
            var pt = new EllPt<FracPoly<FracPoly<ZnInt>>>(X, Y);
            Console.WriteLine(Erl.ContainsPt(pt));
            var pt2 = Erl.Op(pt, pt);
            Console.WriteLine(pt2);
            Console.WriteLine(Erl.ContainsPt(pt2));
            var nPt = Erl.Times(pt, p % l);

            Console.WriteLine(new { nPt });
            Console.WriteLine(Erl.ContainsPt(nPt));
            var phi = FrobRl(pt, 1, psi, psi.Zero, lhs.Num, rhs.Num);
            Console.WriteLine(new { phi });
            Console.WriteLine(Erl.ContainsPt(phi));
            var phi2 = FrobRl(pt, 2, psi, psi.Zero, lhs.Num, rhs.Num);
            Console.WriteLine(new { phi2 });
            Console.WriteLine(Erl.ContainsPt(phi2));
            Console.WriteLine();
            var sum = Erl.Op(phi2, nPt);
            Console.WriteLine(new { sum });
            Console.WriteLine(Erl.ContainsPt(sum));
            Console.WriteLine();
            foreach (var t in l.SeqLazy(1))
            {
                var t_phi = Erl.Times(phi, -t);
                Console.WriteLine(Erl.ContainsPt(t_phi));
                var test = Erl.Op(sum, t_phi);
                Console.WriteLine(new { t, t_phi, test });
                if (test.IsO)
                {
                    frobTr[l] = AmodP(-t, l);
                    break;
                }
            }

            if (!frobTr.ContainsKey(l))
                throw new();
        }

        frobTr.Println("Frob Traces");
        var keys = frobTr.Keys.ToArray();
        var values = keys.Select(k => frobTr[k]).ToArray();
        var crtTable = NumberTheory.CrtTable(keys);
        var L = keys.Aggregate((li, lj) => li * lj);
        var crt = NumberTheory.CRT(values, crtTable, L);
        var ap1 = crt < L / 2 ? crt : crt - L;
        Console.WriteLine($"p = {p} ap = {ap} crt = {crt} ap1 = {ap1} L = {L} Check:{ap == ap1}");
        Console.WriteLine();
    }
}

{
    // EllApFrobTrace([1, 0]);
    RingRl([1, 0]);
}