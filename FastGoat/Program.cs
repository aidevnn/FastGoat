using System.Collections;
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
using RlElt = FastGoat.UserGroup.EllCurve.Frac<FastGoat.UserGroup.EllCurve.Frac<FastGoat.UserGroup.Integers.ZnInt>>;


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

EllPt<GFelt> FrobEq(EllGroup<GFelt> E, EllPt<GFelt> pt, int l, int t)
{
    if (pt.IsO)
        return pt;

    var p = pt.X.P;
    var phi2 = Frob(Frob(pt));
    var t_phi = Frob(E.Times(pt, t));
    var q = E.Times(pt, p % l);
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

int FrobTrace(EllGroup<GFelt> E, EllPt<GFelt>[] nTors, int l) =>
    l.SeqLazy().First(t => nTors.All(pt => Relation(E, pt, t)));

(Dictionary<int, KPoly<Rational>> divPolys,
    Dictionary<int, Dictionary<int, (GFelt g, EllPt<GFelt> e1, EllPt<GFelt> e2)>> allBasis)
    EllApFrobTrace(BigInteger[] curve)
{
    GlobalStopWatch.AddLap();
    var E = EC.EllCoefs(curve);
    var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
    var Ell = E.ToEllGroup();
    var allList = SmallPrimesList(N);
    var pmax = allList.Keys.Max();
    var lmax = allList.Max(e => e.Value.Max());
    Console.WriteLine($"{Ell}");

    var divPolys = EC.DivisionPolynomial(Ell, lmax + 3).f.ToDictionary(e => e.Key, e => e.Value.PrimitiveZPoly());
    divPolys.Println("divPolys");

    Console.WriteLine($"N = {N} pmax = {pmax} listMax = {lmax}");
    var allBasis = new Dictionary<int, Dictionary<int, (GFelt g, EllPt<GFelt> e1, EllPt<GFelt> e2)>>();
    foreach (var (p, listL) in allList)
    {
        if (p < 5)
            continue;

        var ap = EC.EllAp(Ell, p);
        var frobTr = new Dictionary<int, int>();
        var basis = new Dictionary<int, (GFelt g, EllPt<GFelt> e1, EllPt<GFelt> e2)>();
        foreach (var l in listL)
        {
            if (p % l == 0)
                continue;

            var Ep = Ell.ToGF(p);
            Ep.Field = $"GF({p})";
            Console.WriteLine($"{l}-torsion of {Ep} from {Ell}");
            Console.WriteLine($"{l}-divPol mod {p} = {divPolys[l].ToZnPoly(p)}");
            var (nTors, g) = EC.EllFpNTors(Ell, divPolys[l], p, l);
            var Egf = Ell.ToGF(g);
            var (e1, e2) = nTors.GetGenerators().Deconstruct();
            frobTr[l] = FrobTrace(Egf, nTors.GetGenerators().ToArray(), l);
            basis[l] = (g, e1, e2);
            setGF.Add((p, g.X.F.Degree));
        }

        allBasis[p] = basis;
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
    return (divPolys, allBasis);
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
    Console.WriteLine($"{act.Method.Name}");
    Console.WriteLine($"Nb Curves {nbCurves}");
    Console.WriteLine($"Total Errors {nbErrors}/{nbAp}");
    GlobalStopWatch.Show();

    if (setGF.Count != 0)
    {
        var (p, n) = setGF.MaxBy(e => BigInteger.Pow(e.p, e.n));
        Console.WriteLine($"Max field GF({p}^{n})");
        setGF.Where(
                e => !PolynomExt.AllConwayPolys.ContainsKey(e.p) || !PolynomExt.AllConwayPolys[e.p].ContainsKey(e.n))
            .Println(gf => $"GF({gf.p}^{gf.n})", "Missing Fq");
    }

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

void Verif(EllGroupSymb E, EllPt<RlElt> pt, EllPt<GFelt> b1, EllPt<GFelt> b2, EllPt<GFelt> pt1, EllPt<GFelt> pt2,
    string lbl = "pt")
{
    Console.WriteLine($"{E} contains {lbl} => {E.ContainsPt(pt)}");
    var ept1 = E.SubstitutePt(pt, b1);
    var ept2 = E.SubstitutePt(pt, b2);
    Console.WriteLine($"Check on base: {ept1.Equals(pt1)}/{ept2.Equals(pt2)}");
    Console.WriteLine($"b1: {ept1} / {pt1}");
    Console.WriteLine($"b2: {ept2} / {pt2}");
}

void LDivRing(BigInteger[] curve, LogLevel log = LogLevel.Level1)
{
    var allBasis =
        new Dictionary<int, Dictionary<int, (EPoly<ZnInt> g, EllPt<EPoly<ZnInt>> e1, EllPt<EPoly<ZnInt>> e2)>>();

    var E = EC.EllCoefs(curve);
    var El = E.ToLongWeierstrassForm();
    var Ell = El.ToEllGroup();
    if (log != LogLevel.Off)
        allBasis = EllApFrobTrace(El.ArrModel.Select(e => e.Num).ToArray()).allBasis;

    GlobalStopWatch.AddLap();
    var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
    var allList = SmallPrimesList(N);
    var pmax = allList.Keys.Max();
    var lmax = allList.Max(e => e.Value.Max());
    Console.WriteLine($"{E.ToEllGroup()} => {Ell}");

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
        var basis = new Dictionary<int, (EPoly<ZnInt> g, EllPt<EPoly<ZnInt>> e1, EllPt<EPoly<ZnInt>> e2)>();
        if (log != LogLevel.Off)
            basis = allBasis[p];

        foreach (var l in listL)
        {
            GlobalStopWatch.AddLap();
            var psi = divPolys[l].ToZnPoly(p).SubstituteChar('X');
            var facts = IntFactorisation.MusserSFF(psi)
                .SelectMany(e => IntFactorisation.CantorZassenhausAECF(e.g, new(p, a0), p))
                .OrderBy(f => f.Degree).ToArray();

            var Erl = EllGroupSymb.FromEllGroup(E.ToLongWeierstrassForm().ToEllGroup().ToZnInt(p), psi);
            Console.WriteLine($"p={p} l={l} {Erl}");
            Console.WriteLine($"{Erl.EllEq} = 0");
            Console.WriteLine($"psi = {psi}");
            Console.WriteLine($"psi ~ {facts.Glue(" * ", "({0})")}");

            var pt = Erl.Pt;
            var p_Pt = Erl.Times(pt, p % l);
            var phi = Erl.FrobRl(pt);
            var phi2 = Erl.FrobRl(pt, 2);
            var add_phi2_p = Erl.Op(phi2, p_Pt);

            if (log != LogLevel.Off)
            {
                var (bg, be1, be2) = basis[l];
                var Ep = Ell.ToGF(bg);

                Console.WriteLine(new { pt });
                Verif(Erl, pt, be1, be2, be1, be2, "P");

                var pti = Erl.Invert(pt);
                Console.WriteLine(new { pti });
                Verif(Erl, pti, be1, be2, Ep.Invert(be1), Ep.Invert(be2), "-P");

                var pt2 = Erl.Op(pt, pt);
                Console.WriteLine(new { pt2 });
                Verif(Erl, pt2, be1, be2, Ep.Times(be1, 2), Ep.Times(be2, 2), "2P");

                Console.WriteLine(new { p_Pt });
                Verif(Erl, p_Pt, be1, be2, Ep.Times(be1, p % l), Ep.Times(be2, p % l), "[p]*P");

                Console.WriteLine(new { phi });
                Verif(Erl, phi, be1, be2, Frob(be1), Frob(be2), "phi");

                Console.WriteLine(new { phi2 });
                Verif(Erl, phi2, be1, be2, Frob(be1, 2), Frob(be2, 2), "phi^2");

                Console.WriteLine(new { add_phi2_pPy = add_phi2_p });
                Verif(Erl, add_phi2_p, be1, be2, Ep.Op(Frob(be1, 2), Ep.Times(be1, p % l)),
                    Ep.Op(Frob(be2, 2), Ep.Times(be2, p % l)), "phi^2 + [p]*P");
            }

            foreach (var t in l.SeqLazy())
            {
                var t_phi = Erl.FrobRl(Erl.Times(pt, -t));
                var eqFrob = Erl.Op(add_phi2_p, t_phi);
                if (log != LogLevel.Off)
                {
                    var (bg, be1, be2) = basis[l];
                    var Ep = Ell.ToGF(bg);
                    Verif(Erl, t_phi, be1, be2, Frob(Ep.Times(be1, -t)), Frob(Ep.Times(be2, -t)), "-[t]*phi");
                    Verif(Erl, eqFrob, be1, be2, FrobEq(Ep, be1, l, -t), FrobEq(Ep, be2, l, -t),
                        "phi^2 - [t]*phi + [p]*P");
                    Console.WriteLine(new { t, eqFrob });
                }

                if (eqFrob.IsO)
                {
                    frobTr[l] = t;
                    Console.WriteLine($"t = {t}");
                    break;
                }
            }

            if (!frobTr.ContainsKey(l))
                throw new();

            GlobalStopWatch.Show();
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
        if (ap != ap1)
        {
            ++nbErrors;
            Console.WriteLine($"Warning[{nbErrors}]");
        }

        ++nbAp;
        Console.WriteLine();
    }

    ++nbCurves;
    GlobalStopWatch.Show($"End LDivRing {Ell}");
    Console.WriteLine();
}

void testLDivRing1()
{
    LDivRing([1, 0]);
    LDivRing([-1, 0]);
    // LDivRing([1, 1]);
    // LDivRing([-43, 166]);
}

void testLDivRing2()
{
    LDivRing([0, 0, 1, 0, -7]);
    LDivRing([1, 0, 0, 0, 1]);
    LDivRing([1, -1, 0, -4, 4]);
    LDivRing([1, -1, 1, -19353, 958713]);
    LDivRing([1, 1, 1, -17714, 900047]);
}

void testLDiv()
{
    runTest(testLDivRing1);
    runTest(testLDivRing2);

    // testLDivRing1
    // Nb Curves 4
    // Total Errors 0/20
    // #  Time:2m26s
    // 
    // testLDivRing2
    // Nb Curves 5
    // Total Errors 0/49
    // #  Time:14m24s
    // 

    // Console.Beep();
}

Polynomial<K, Xi> Eq<K>(Polynomial<K, Xi> x, Polynomial<K, Xi> y, Polynomial<K, Xi> a1, Polynomial<K, Xi> a2,
    Polynomial<K, Xi> a3, Polynomial<K, Xi> a4, Polynomial<K, Xi> a6)
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    => (x.Pow(3) + a2 * x.Pow(2) + a4 * x + a6) - (y * y + a1 * x * y + a3 * y);

void SymbOps()
{
    Ring.DisplayPolynomial = MonomDisplay.StarPowFct;
    var xis = Ring.Polynomial(Rational.KOne(), MonomOrder.Lex, "x1", "x2", "x3", "y1", "y2", "y3", "a1", "a2", "a3",
        "a4", "a6", "alpha");
    var (x1, x2, x3, y1, y2, y3) = xis.Deconstruct();
    var (a1, a2, a3, a4, a6, alpha) = xis.Skip(6).Deconstruct();

    var eq1 = Eq(x1, y1, a1, a2, a3, a4, a6);
    var eq3 = Eq(x3, y3, a1, a2, a3, a4, a6);
    var eq2 = eq3.Substitute(alpha * (x3 - x1) + y1, "y3");
    var dec = Ring.Decompose(eq2, "x3").Item1;
    Console.WriteLine($"{eq3} = 0");
    dec.Println($"{eq2} = 0");
    Console.WriteLine();

    Console.WriteLine("case -P1 = P3");
    Console.WriteLine($"var x3 = {x1};");
    Console.WriteLine($"var y3 = {-a1 * x1 - a3 - y1};");
    Console.WriteLine();

    var sol = -(dec[x3 * x3] + x1 + x2);
    Console.WriteLine("case P1 + P2 = P3");
    Console.WriteLine($"var {alpha} = ({y2 - y1}) / ({x2 - x1});");
    Console.WriteLine($"var x3 = {sol};");
    Console.WriteLine($"var y3 = {alpha * (x3 - x1) + y1};");
    Console.WriteLine();

    Console.WriteLine("case 2*P1 = P3");
    Console.WriteLine($"var {alpha} = ({eq1.D("x1")}) / ({-eq1.D("y1")});");
    Console.WriteLine($"var x3 = {sol.Substitute(x1, "x2")};");
    Console.WriteLine($"var y3 = {alpha * (x3 - x1) + y1};");
    Console.WriteLine();
}

void testOps1()
{
    var E = EC.EllGroup([7, 0, 0, 16, 0]);
    Console.WriteLine(E.Eq);
    var O = E.O;
    EllPt<Rational> P = ("-8", "40");
    EllPt<Rational> Q = ("4", "4");
    Console.WriteLine(new { O, P, Q });
    Console.WriteLine($"-P = {E.Invert(P)}");
    Console.WriteLine($"-Q = {E.Invert(Q)}");
    Console.WriteLine($"P + Q = {E.Op(P, Q)}");
    Console.WriteLine($"2P = {E.Op(P, P)}");
    Console.WriteLine($"2P = {E.Times(P, 2)}");
    Console.WriteLine($"2Q = {E.Times(Q, 2)}");
    Console.WriteLine($"2P + 2Q = {E.Op(E.Times(P, 2), E.Times(Q, 2))}");
    Console.WriteLine($"2(P + Q) = {E.Times(E.Op(P, Q), 2)}");
    Console.WriteLine();
}

void testOps2()
{
    var E = EC.EllGroup([1, -1, 1, -14, 29]);
    Console.WriteLine(E.Eq);
    var O = E.O;
    EllPt<Rational> P = ("-3", "7");
    EllPt<Rational> Q = ("9", "19");
    Console.WriteLine(new { O, P, Q });
    Console.WriteLine($"-P = {E.Invert(P)}");
    Console.WriteLine($"-Q = {E.Invert(Q)}");
    Console.WriteLine($"P + Q = {E.Op(P, Q)}");
    Console.WriteLine($"2P = {E.Op(P, P)}");
    Console.WriteLine($"2P = {E.Times(P, 2)}");
    Console.WriteLine($"2Q = {E.Times(Q, 2)}");
    Console.WriteLine($"2P + 2Q = {E.Op(E.Times(P, 2), E.Times(Q, 2))}");
    Console.WriteLine($"2(P + Q) = {E.Times(E.Op(P, Q), 2)}");
    Console.WriteLine();
}

ConcreteGroup<EllPt<GFelt>> EllFq(BigInteger[] curve, int q)
{
    var g = NumberTheory.PrimitiveRoot(FG.FqX(q, 'a'));
    EllGroup<GFelt> E = EC.EllGroup(curve).ToGF(g);
    if (q > 3000)
        throw new();

    var d = 2 * double.Sqrt(q);
    var y = FG.KPoly('y', g);
    var gens = new HashSet<EllPt<GFelt>>();
    var set = new HashSet<EllPt<GFelt>>();
    set.Add(E.Neutral());
    var idx = 0;
    foreach (var x in q.SeqLazy().Shuffle().Select(k => g.Pow(k)).Prepend(g.Zero))
    {
        if (set.Any(pt => !pt.IsO && pt.X.Equals(x)))
            continue;

        var eq = (y.Pow(2) + E.a1 * x * y + E.a3 * y) - (x.Pow(3) + E.a2 * x.Pow(2) + E.a4 * x + E.a6);
        var pts = IntFactorisation.FirrFsepBerlekampAECF(eq, g, q)
            .Where(f => f.g.Degree == 1)
            .Select(f => (x, y: -f.g[0]))
            .Select(e => new EllPt<GFelt>(e.x, e.y))
            .Where(e => E.Contains(e.X, e.Y))
            .ToList();

        idx++;
        if (pts.Count == 0)
            continue;

        gens.UnionWith(pts);
        set = Group.GenerateElements(E, set, pts);
        var count = set.Count;
        // Console.WriteLine(new { x, newGens = pts.Count, ptsCount = count });
        if (idx > g.P && idx > d && double.Abs(count - q - 1) < d)
            break;
    }

    return Group.Generate(E, gens.ToArray());
}

// Mathematics of Public Key Cryptography. Version 2.0
// Steven D Galbraith
// CHAPTER 9. ELLIPTIC CURVES
// 9.10. FROBENIUS MAP
// Corollary 9.10.9 page 185
(EPoly<Rational> alpha, EPoly<Rational> beta) SolveFrobTrace(EllGroup<Rational> E0, int p)
{
    var t = EC.EllAp(E0, p);
    var x = FG.QPoly();
    var (alpha, beta) = IntFactorisation.AlgebraicRoots(x.Pow(2) - t * x + p).Deconstruct();
    return (alpha, beta);
}

int EllCardExt(EPoly<Rational> alpha, EPoly<Rational> beta, int q)
{
    var p = Primes10000.First(p => q % p == 0);
    var n = (int)double.Round(double.Log(q) / double.Log(p));
    var card = (1 - alpha.Pow(n)) * (1 - beta.Pow(n));
    return (int)card[0].Num;
}

void testEllFq(BigInteger[] curve, int p)
{
    var n = (int)(double.Log(1000) / double.Log(p));
    var E = EC.EllGroup(curve);
    var (alpha, beta) = SolveFrobTrace(E, p);
    for (int i = 1; i <= n; i++)
    {
        var q = p.Pow(i);
        var card = EllCardExt(alpha, beta, q);
        Console.WriteLine($"q={q,-6} #{$"{E.ToGF(q)}",-25} = {card}");
        var gEll = EllFq(curve, q);
        var abEll = Group.AbelianDecompositions(gEll);
        var abType = abEll.abType.Select(e => e.o).Glue(" x ");
        abEll.abType.Println(l => $"{l.g} of order {l.o}", $"Generators of {gEll.ShortName} ~ [{abType}]");
        if (card != gEll.Count())
            throw new();

        Console.WriteLine();
    }
}

void testEllFqCard(BigInteger[] curve, int p, bool bf = false)
{
    var n = (int)(double.Log(100000) / double.Log(p));
    var E = EC.EllGroup(curve);
    var (alpha, beta) = SolveFrobTrace(E, p);
    for (int i = 1; i <= n; i++)
    {
        var q = p.Pow(i);
        var cardExt = EllCardExt(alpha, beta, q);
        if (q > 3000 || !bf)
            Console.WriteLine($"q={q,-6} #{$"{E.ToGF(q)}",-25} = {cardExt}");
        else
        {
            var cardBf = EllCardBF(E, q);
            Console.WriteLine($"q={q,-6} #{$"{E.ToGF(q)}",-25} = {cardExt,-7} cardBf = {cardBf}");
            if (cardExt != cardBf)
                throw new();
        }
    }

    Console.WriteLine();
}

IEnumerable<GFelt> SetFq(GFelt g)
{
    var (p, n) = (g.P, g.F.Degree);
    yield return g.Zero;
    var x = g.One;
    for (int i = 0; i < p.Pow(n) - 1; i++)
    {
        yield return x;
        x *= g;
    }
}

int EllCardBF(EllGroup<Rational> E0, int q)
{
    var card = 1;
    var g = NumberTheory.PrimitiveRoot(FG.FqX(q, 'a'));
    EllGroup<GFelt> E = E0.ToGF(g);
    if (int.IsPow2(q))
    {
        // Mathematics of Public Key Cryptography. Version 2.0
        // Steven D Galbraith
        // CHAPTER 2. BASIC ALGORITHMIC NUMBER THEORY
        // 2.14.2 Solving Quadratic Equations in Finite Fields
        // Exercise 2.14.7. page 67
        foreach (var x in SetFq(g))
        {
            var a = E.a1 * x + E.a3;
            var b = x.Pow(3) + E.a2 * x.Pow(2) + E.a4 * x + E.a6;
            if (a.IsZero())
            {
                if (b.IsZero() || NumberTheory.LegendreJacobiGf(b).IsOne())
                    ++card;
            }
            else
            {
                if (IntFactorisation.Mk(b / a.Pow(2), 1, q).IsZero())
                    card += 2;
            }
        }
    }
    else
    {
        var Eql = E.ToLongWeierstrassForm();
        foreach (var x in SetFq(g))
        {
            var resQuad = NumberTheory.LegendreJacobiGf(x.Pow(3) + Eql.a2 * x.Pow(2) + Eql.a4 * x + Eql.a6);
            if (resQuad.IsOne())
                card += 2;
            else if (resQuad.IsZero())
                ++card;
        }
    }

    return card;
}

{
    // testLDiv();
    // SymbOps();
    // testOps1();
    // testOps2();

    BigInteger[][] curves =
    [
        [0, 1, 1, 0, 1],
        [0, 0, 1, 0, 1],
        [1, 0, 1, 0, 1],
        [1, 1, 0, 0, 1],
        [1, 0, 0, 1, 0],
        [1, 0, 2, 0, 1],
        [2, 0, 0, 1, 0]
    ];

    foreach (var p in new[] { 2, 3, 5 })
    {
        foreach (var curve in curves.Where(e => !EC.EllCoefs(e).ToZnInt(p).Disc.IsZero()))
            testEllFq(curve, p);
    }

    foreach (var p in Primes10000.Take(10))
    {
        foreach (var curve in curves.Where(e => !EC.EllCoefs(e).ToZnInt(p).Disc.IsZero()))
            testEllFqCard(curve, p, bf: true);
    }
}