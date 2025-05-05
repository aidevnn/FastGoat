using System.Collections;
using System.Numerics;
using System.Reflection.Metadata.Ecma335;
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
using RlElt = FastGoat.Structures.VecSpace.Frac<FastGoat.Structures.VecSpace.Frac<FastGoat.UserGroup.Integers.ZnInt>>;


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

Dictionary<int, Dictionary<int, (GFelt g, EllPt<GFelt> e1, EllPt<GFelt> e2)>> EllApFrobTrace(BigInteger[] curve)
{
    GlobalStopWatch.AddLap();
    var E = EC.EllCoefs(curve);
    var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
    var Ell = E.ToEllGroup();
    var allList = SmallPrimesList(N);
    var pmax = allList.Keys.Max();
    var lmax = allList.Max(e => e.Value.Max());
    Console.WriteLine($"{Ell} {Ell.ShortFormStr} {Ell.LongFormStr}");

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
    return allBasis;
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
        setGF.Where(e => !PolynomExt.AllConwayPolys.ContainsKey(e.p) || !PolynomExt.AllConwayPolys[e.p].ContainsKey(e.n))
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

void Verif(EllGroupSymb E,  EllPt<RlElt> pt, EllPt<GFelt> b1, EllPt<GFelt> b2, EllPt<GFelt> pt1, EllPt<GFelt> pt2,
    string lbl = "pt")
{
    Console.WriteLine($"{E} contains {lbl} => {E.ContainsPt(pt)}");
    var ept1 = E.SubstitutePt(pt, b1);
    var ept2 = E.SubstitutePt(pt, b2);
    Console.WriteLine($"Check on base: {ept1.Equals(pt1)}/{ept1.Equals(pt2)}");
    Console.WriteLine($"b1: {ept1} / {pt1}");
    Console.WriteLine($"b2: {ept2} / {pt2}");
}

void LDivRing(BigInteger[] curve, LogLevel log = LogLevel.Off)
{
    var allBasis =
        new Dictionary<int, Dictionary<int, (EPoly<ZnInt> g, EllPt<EPoly<ZnInt>> e1, EllPt<EPoly<ZnInt>> e2)>>();

    var E = EC.EllCoefs(curve);
    var El = E.ToLongWeierstrassForm();
    var Ell = El.ToEllGroup();
    if (log != LogLevel.Off)
        allBasis = EllApFrobTrace(El.ArrModel.Select(e => e.Num).ToArray());

    GlobalStopWatch.AddLap();
    var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
    var allList = SmallPrimesList(N);
    var pmax = allList.Keys.Max();
    var lmax = allList.Max(e => e.Value.Max());
    Console.WriteLine($"{E.ToEllGroup()} => {Ell} => {Ell.ShortFormStr}");

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
    LDivRing([1, 1]);
    LDivRing([-43, 166]);
}

void testLDivRing2()
{
    LDivRing([0, 0, 1, 0, -7]);
    LDivRing([1, 0, 0, 0, 1]);
    LDivRing([1, -1, 0, -4, 4]);
    LDivRing([1, -1, 1, -19353, 958713]);
    LDivRing([1, 1, 1, -17714, 900047]);
}

{
    runTest(testLDivRing1);
    runTest(testLDivRing2);
    
    // testLDivRing1
    // Nb Curves 4
    // Total Errors 0/20
    // #  Time:8m42s
    // 
    // testLDivRing2
    // Nb Curves 5
    // Total Errors 0/49
    // #  Time:15m48s
    // 

    Console.Beep();
}

public struct EllGroupSymb : IGroup<EllPt<RlElt>>
{
    public EllGroupSymb(EllCoefs<RlElt> ellCoefs, KPoly<ZnInt> divPol)
    {
        Coefs = ellCoefs.Model;
        Disc = ellCoefs.Disc;

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
        (X, Y) = FG.BiVariateFracZnInt(Disc.P, 'X', 'Y');
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

        var P1 = ConvertToShort(P);
        return ConvertFromShort(new EllPt<RlElt>(P1.X, -P1.Y));
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

        var (e1_, e2_) = (ConvertToShort(e1), ConvertToShort(e2));
        var (x1, y1, x2, y2) = (e1_.X, e1_.Y, e2_.X, e2_.Y);
        if (!SimplifyEll(x1 - x2, DivPolynomial, EllEq.Num).IsZero())
        {
            var alpha = (y2 - y1) / (x2 - x1);
            var x3 = alpha.Pow(2) - x1 - x2;
            var y3 = -y1 + alpha * (x1 - x3);
            return ConvertFromShort(new(x3, y3));
        }
        else
        {
            if (!SimplifyEll(y1 + y2, DivPolynomial, EllEq.Num).IsZero())
            {
                if (y1.IsZero())
                {
                    Console.WriteLine(new { y1, y2 });
                    return new();
                }

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
        var Y = FG.BiVariateFracZnInt(divPoly.P, 'X', 'Y').Y;
        var (a1, a2, a3, a4, a6) = E.Coefs;
        var (A1, A2, A3, A4, A6) = (a1 * Y.KOne * Y.One, a2 * Y.KOne * Y.One, a3 * Y.KOne * Y.One, a4 * Y.KOne * Y.One,
            a6 * Y.KOne * Y.One);
        return new(A1, A2, A3, A4, A6, divPoly);
    }
}