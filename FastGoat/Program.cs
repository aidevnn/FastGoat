using System.Numerics;
using System.Text.Json.Serialization;
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

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

(KPoly<Rational> X2, KPoly<Rational> Y2) EllP2(EllCoefs<Rational> E)
{
    var (a1, a2, a3, a4, a6) = E.Model;
    var xis = Ring.Polynomial(Rational.KZero(), MonomOrder.Lex, "Y", "X");
    var (Y0, X0) = xis.Deconstruct();
    var basis = new PolynomialBasis<Rational, Xi>((Y0.Pow(2) + a1 * X0 * Y0 + a3 * Y0) -
                                                  (X0.Pow(3) + a2 * X0.Pow(2) + a4 * X0 + a6));

    var (Y, X) = Ring.EPolynomial(xis).Deconstruct();
    var Ell = new EllGroup<EPolynomial<Rational>>(a1 * X.One, a2 * X.One, a3 * X.One, a4 * X.One, a6 * X.One);
    var P = new EllPt<EPolynomial<Rational>>(X, Y);

    var P1 = Ell.ConvertToShort(P);
    var (x1, y1) = P1;
    var alpha = (3 * x1.Pow(2) + Ell.ShortForm.A) / (2 * y1);
    var x2 = alpha.Pow(2) - 2 * x1;
    var y2 = -y1 + alpha * (x1 - x2);
    var P2 = Ell.ConvertFromShort(new(x2, y2));
    var xn = basis.Rem(P2.X.Num);
    var yn = basis.Rem(P2.Y.Num);
    P2 = new(new(xn, P2.X.Denom), new(yn.Monic(), P2.Y.Denom / yn.LeadingDetails.lc));
    Console.WriteLine($"{E.Eq}, P = {P} and 2P = {P2}");
    return (xn.ToKPoly("X").Monic, yn.ToKPoly("Y").Monic);
}

KPoly<KPoly<Rational>> Simplify(KPoly<KPoly<Rational>> P, KPoly<KPoly<Rational>> R)
{
    if (P.Degree < 2)
        return P;

    var (quo, rem) = P.Div(P.X.Pow(2));
    return Simplify(quo * R + rem, R);
}

(KPoly<KPoly<Rational>> R, Dictionary<int, KPoly<KPoly<Rational>>> psi, Dictionary<int, KPoly<Rational>> f)
    DivisionPolynomial(EllGroup<Rational> E, int nmax)
{
    var x = FG.QPoly('X');
    var Y = FG.KPoly('Y', x);
    var X = x * Y.One;
    var (A, B) = (E.ShortForm.A * x.One * Y.One, E.ShortForm.B * x.One * Y.One);
    var R = X.Pow(3) + A * X + B;

    var psi = new Dictionary<int, KPoly<KPoly<Rational>>>();
    (psi[0], psi[1], psi[2]) = (Y.Zero, Y.One, 2 * Y);
    psi[3] = 3 * X.Pow(4) + 6 * A * X.Pow(2) + 12 * B * X - A.Pow(2);
    psi[4] = 4 * Y * (X.Pow(6) + 5 * A * X.Pow(4) + 20 * B * X.Pow(3) - 5 * A.Pow(2) * X.Pow(2) - 4 * A * B * X -
                      8 * B.Pow(2) - A.Pow(3));

    for (int n = 2; n < nmax / 2; n++)
    {
        if (n > 2)
            psi[2 * n] = Simplify(psi[n] * (psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2)) / (2 * Y),
                R);
        psi[2 * n + 1] = Simplify(psi[n + 2] * psi[n].Pow(3) - psi[n + 1].Pow(3) * psi[n - 1], R);
    }

    var f = psi.ToDictionary(e => e.Key, e => e.Key % 2 == 0 ? (e.Value / Y)[0].Clone() : e.Value[0].Clone());
    return (R, psi, f);
}

(EllPt<FracPoly<FracPoly<Rational>>> P, EllPt<FracPoly<FracPoly<Rational>>> nP)
    NPt(int n, Dictionary<int, KPoly<KPoly<Rational>>> psi, KPoly<KPoly<Rational>> R)
{
    var Y = FG.KFracPoly('Y', FG.KFracPoly(R.KOne.X).X).X;
    var X = Y.KOne.X + Y.Zero;

    var tmp1 = Simplify(psi[n - 1] * psi[n + 1], R).ToFrac(Y) / Simplify(psi[n].Pow(2), R).ToFrac(Y);
    var nPX = X - tmp1;

    var denom = n % 2 == 0
        ? Y * Simplify(4 * psi[n].Pow(3), R).ToFrac(Y)
        : Simplify(4 * R.X * psi[n].Pow(3), R).ToFrac(Y);
    var tmp2 = Simplify(psi[n + 2] * psi[n - 1].Pow(2) - psi[n - 2] * psi[n + 1].Pow(2), R).ToFrac(Y);
    var nPY = tmp2 / denom;

    return (new EllPt<FracPoly<FracPoly<Rational>>>(X, Y), new EllPt<FracPoly<FracPoly<Rational>>>(nPX, nPY));
}

EllPt<BigCplx>[] EllCM(EllGroup<Rational> E, KPoly<Rational> EqX)
{
    var (A, B) = (E.ShortForm.A, E.ShortForm.B);
    var solsX = FG.NRoots(EqX.ToBcPoly(60), maxLoop: 800).ToHashSet();
    var pts = new HashSet<EllPt<BigCplx>>();
    foreach (var x in solsX)
    {
        var y2 = x.Pow(3) + A * x + B;
        var y = BigCplx.Sqrt(y2).ToBigCplx(10);
        pts.Add(new(x.ToBigCplx(10), y));
        pts.Add(new(x.ToBigCplx(10), -y));
    }

    return pts.ToArray();
}

EPoly<ZnInt>[] Roots(KPoly<EPoly<ZnInt>> P, EPoly<ZnInt> a, BigInteger q)
{
    var facts = IntFactorisation.MusserSFF(P)
        .SelectMany(f => IntFactorisation.CantorZassenhausAECF(f.g, a, q))
        .ToArray();
    if (facts.All(f => f.Degree == 1))
        return facts.Select(f => -f[0]).ToArray();

    return [];
}

ConcreteGroup<EllPt<EPoly<ZnInt>>>
    EpNTors(EllGroup<Rational> E, KPoly<Rational> P, int p, int l, LogLevel log = LogLevel.Level1)
{
    if (log != LogLevel.Off)
        GlobalStopWatch.AddLap();

    EllGroup<EPoly<ZnInt>> Efq = E.ToGF(p);
    var setPts = new HashSet<EllPt<EPoly<ZnInt>>>();

    for (int i = 1; i < 10; i++)
    {
        var q = BigInteger.Pow(p, i);

        Efq = E.ToGF(q);
        var Pq = P.ToGF(q);
        Efq.Field = $"GF({p}^{i})";
        var a = Pq.KOne.X;
        var tmpX = Roots(Pq, a, q);
        if (tmpX.Length == 0)
        {
            if (log == LogLevel.Level2)
                Console.WriteLine($"P = {Pq} dont split in GF({p}^{i})[X]");
            continue;
        }

        if (log == LogLevel.Level2)
        {
            Console.WriteLine($"P = {Pq} split in GF({p}^{i})[X]");
            Console.WriteLine($"x=[{tmpX.Glue(", ")}]");
        }

        var (A, B) = (Efq.ShortForm.A, Efq.ShortForm.B);
        var X = FG.KPoly('X', Efq.a3.X);
        var Y = FG.KPoly('Y', Efq.a3.X);
        var tmpY = new Dictionary<EPoly<ZnInt>, EPoly<ZnInt>[]>();
        foreach (var x in tmpX)
        {
            var sols = Roots(Y.Pow(2) - (x.Pow(3) + A * x + B), a, q);
            if (sols.Length > 0)
            {
                if (log == LogLevel.Level2)
                    Console.WriteLine($"x={x} y=[{sols.Glue(", ")}]");
                tmpY[x] = sols;
            }
            else
            {
                tmpY.Clear();
                break;
            }
        }

        if (tmpY.Count == 0)
        {
            if (log == LogLevel.Level2)
                Console.WriteLine($"Y^2 = {X.Pow(3) + A * X + B} dont have {l}-torsion in GF({p}^{i})");
            continue;
        }

        setPts.UnionWith(tmpY.SelectMany(e => e.Value.Select(y => new EllPt<EPoly<ZnInt>>(e.Key, y))));

        if (log == LogLevel.Level2)
            setPts.Println("Points");

        break;
    }

    var gEll = Group.Generate(Efq, setPts.Where(e => Efq.Contains(e.X, e.Y)).ToArray());
    var Cn = Group.Generate($"{l}-Tors({gEll})", gEll,
        gEll.ElementsOrders.Where(e => e.Value == l).Select(e => e.Key).ToArray());

    if (log == LogLevel.Level2)
        DisplayGroup.HeadElements(gEll);

    if (log != LogLevel.Off)
    {
        DisplayGroup.HeadElements(Cn);

        var abCn = FG.AbelianDirectSum(Cn);
        abCn.DecompMap.Println(e => $"{e.Key} of order {e.Value}", $"Generators of {Cn}");
        Console.WriteLine();
        Console.WriteLine($"{Cn} ~ {abCn.DecompMap.Values.Glue(" x ", "C{0}")}");

        if (Cn.Count() == 1)
            Console.WriteLine($"Warnings");

        GlobalStopWatch.Show();
        Console.WriteLine();
    }

    return Cn;
}

void CplxMul_NTors(BigInteger[] curve)
{
    var E = EC.EllCoefs(curve);
    Console.WriteLine($"Ell({E.ModelStr})(Q) {E.Eq}");
    var (R, psi, divPolys) = DivisionPolynomial(E.ToEllGroup(), 8);
    var (Ell, N, _) = EC.EllTateAlgorithm(E);
    var nmax = 2 * double.Sqrt(N);

    EllP2(E);
    var (P, P2) = NPt(2, psi, R);
    var P3 = NPt(3, psi, R).nP;
    var P4 = NPt(4, psi, R).nP;
    var P5 = NPt(5, psi, R).nP;
    Console.WriteLine($"P = {P}");

    Console.WriteLine($"2P = {P2}");
    Console.WriteLine(P2.X);
    var x2P = P2.X.Num[0].Denom.Monic;
    Console.WriteLine(x2P);
    var facts2P = IntFactorisation.FactorsQ(x2P);
    facts2P.Println($"factors x(2P) = {x2P}");
    // EllCM(Ell, x2P).Println("C[2]");
    Console.WriteLine();

    Console.WriteLine($"3P = {P3}");
    Console.WriteLine(P2.X - P.X);
    Console.WriteLine(P3.X);
    var x3Pa = P3.X.Num[0].Denom.Monic;
    var x3Pb = (P2.X - P.X).Num[0].Num.Monic;
    Console.WriteLine(x3Pa);
    Console.WriteLine(x3Pb);
    Console.WriteLine(x3Pa.Div(x3Pb));
    var facts3P = IntFactorisation.FactorsQ(x3Pb);
    facts3P.Println($"factors x(3P) = {x3Pb}");
    // EllCM(Ell, x3Pa).Println("C[3]");
    Console.WriteLine();

    Console.WriteLine($"4P = {P4}");
    Console.WriteLine(P3.X - P.X);
    Console.WriteLine(P4.X);
    var x4Pa = P4.X.Num[0].Denom.Monic;
    var x4Pb = (P3.X - P.X).Num[0].Num.Monic;
    var x4Pc = (P2.Y * P.X).Num[0].Num.Monic;
    Console.WriteLine(x4Pa);
    Console.WriteLine(x4Pb);
    Console.WriteLine(x4Pc);
    Console.WriteLine();
    Console.WriteLine(x4Pa.Div(x4Pc));
    Console.WriteLine(x4Pb.Div(x4Pc));
    var facts4P = IntFactorisation.FactorsQ(x4Pc);
    facts4P.Println($"factors x(4P) = {x4Pc}");
    // EllCM(Ell, x4Pa).Println("C[4]");
    Console.WriteLine();

    var x5Pa = P5.X.Num[0].Denom.Monic;
    Console.WriteLine(x5Pa);
    var facts5P = IntFactorisation.FactorsQ(x5Pa);
    facts5P.Println($"factors x(5P) = {x5Pa}");
    var x5Pb = facts5P.Aggregate(x5Pa.One, (acc, fi) => acc * fi.Item1);
    Console.WriteLine(x5Pb);
    Console.WriteLine();

    foreach (var p in Primes10000.Where(p => p > 3 && p <= nmax && !E.Disc.Mod(p).IsZero()))
    {
        EpNTors(Ell, x2P, p, 2);
        EpNTors(Ell, x3Pb, p, 3);
        EpNTors(Ell, x4Pc, p, 4);
        if (p > 5)
            EpNTors(Ell, x5Pb, p, 5);
    }
}

(List<int> listL, int L) GetPrimesL(int p)
{
    var listL = new List<int>();
    var L = 1;
    var lmax = 8 * double.Sqrt(p);
    foreach (var l in Primes10000.Where(l => l < p))
    {
        listL.Add(l);
        L *= l;
        if (L > lmax)
            break;
    }

    return (listL, L);
}

(string ModelStr, string N, int pmax, int[] listLmax, int lmax) ParamsNTors(BigInteger[] curve, Rational N)
{
    var E = EC.EllCoefs(curve);
    var nmax = 2 * double.Sqrt(N);
    var pmax = Primes10000.Last(p => p <= nmax && !N.Mod(p).IsZero());
    var (listLmax, Lmax) = GetPrimesL(pmax);

    Console.WriteLine($"Ell({E.ModelStr})(Q) pmax = {pmax} lmax = {listLmax.Max()}  L = {Lmax}");
    return (E.ModelStr, N: $"{N}", pmax, listLmax.ToArray(), lmax: listLmax.Max());
}

EllPt<EPoly<ZnInt>> Frob(EllPt<EPoly<ZnInt>> pt)
{
    if (pt.IsO)
        return pt;

    var p = pt.X.P;
    return new(pt.X.Pow(p), pt.Y.Pow(p));
}

bool Relation(ConcreteGroup<EllPt<EPoly<ZnInt>>> E, EllPt<EPoly<ZnInt>> pt, int t)
{
    if (pt.IsO)
        return true;

    var p = pt.X.P;
    var phi2 = Frob(Frob(pt));
    var t_phi = Frob(E.Times(pt, -t));
    var q = E.Times(pt, p);
    return E.Op(phi2, E.Op(t_phi, q)).IsO;
}

int FrobTrace(ConcreteGroup<EllPt<EPoly<ZnInt>>> nTors, int l)
{
    var p = nTors.First(pt => !pt.IsO).X.P;
    return (p - 1).SeqLazy(1).First(t => nTors.All(pt => Relation(nTors, pt, t))) % l;
}

void EllApFrobTrace(BigInteger[] curve, Rational N)
{
    var E = EC.EllCoefs(curve);
    var Ell = E.ToEllGroup();
    var nmax = 2 * double.Sqrt(N);

    var lmax = ParamsNTors(curve, N).lmax;
    var (R, psi, divPolys) = DivisionPolynomial(Ell, 2 * lmax);
    var prepXnP = Primes10000.Where(l => l <= lmax).ToDictionary(l => l, l => NPt(l, psi, R).nP.X.Num[0].Denom.Monic);

    foreach (var p in Primes10000.Where(p => p > 3 && p <= nmax && !N.Mod(p).IsZero()))
    {
        var Ep = Ell.ToGF(p);
        Ep.Field = $"GF({p})";
        var (listL, _) = GetPrimesL(p);
        var ap = EC.EllAp(Ell, p);
        var frobTr = new Dictionary<int, int>();
        foreach (var l in listL)
        {
            if (p % l == 0)
                continue;

            var nTors = EpNTors(Ell, prepXnP[l], p, l);
            if (nTors.Count() != l * l)
                continue;

            frobTr[l] = FrobTrace(nTors, l);
            Console.WriteLine();
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

void testDivPolys()
{
    var E = EC.EllCoefs([1, 0]);
    Console.WriteLine($"Ell({E.ModelStr})(Q) {E.Eq}");
    var (R, psi, divPolys) = DivisionPolynomial(E.ToEllGroup(), 8);
    psi.Println("psi");
    divPolys.Println("divPolys");
    divPolys.Where(e => e.Key != 0)
        .ToDictionary(e => e.Key, e => (e.Value.Degree, (e.Key.Pow(2) - 1 - 3 * ((e.Key + 1) % 2)) / 2))
        .Println("Degree of divPolys");

    EllP2(E);
    var (P, P2) = NPt(2, psi, R);
    var P3 = NPt(3, psi, R).nP;
    Console.WriteLine($"P = {P}");
    Console.WriteLine($"2P = {P2}");
    Console.WriteLine($"3P = {P3}");
}

void runCM_NTors()
{
    GlobalStopWatch.Restart();

    GlobalStopWatch.AddLap();
    CplxMul_NTors([1, 0]);
    GlobalStopWatch.Show();
    Console.WriteLine();

    GlobalStopWatch.AddLap();
    CplxMul_NTors([-1, 0]);
    GlobalStopWatch.Show();
    Console.WriteLine();

    GlobalStopWatch.AddLap();
    CplxMul_NTors([1, 1]);
    GlobalStopWatch.Show();
    Console.WriteLine();

    Console.Beep();
}

{
    GlobalStopWatch.Restart();

    {
        BigInteger[] curve = [1, 0];
        var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
        EllApFrobTrace(curve, N);
    }
    {
        BigInteger[] curve = [-1, 0];
        var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
        EllApFrobTrace(curve, N);
    }
    {
        BigInteger[] curve = [1, 1];
        var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
        EllApFrobTrace(curve, N);
    }

    GlobalStopWatch.Show();
}