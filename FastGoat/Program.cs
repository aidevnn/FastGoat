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

EllPt<GFelt> FrobEq(EllGroup<GFelt> E, EllPt<GFelt> pt, int l, int t)
{
    if (pt.IsO)
        return pt;

    var p = pt.X.P;
    var phi2 = EC.Frob(EC.Frob(pt));
    var t_phi = EC.Frob(E.Times(pt, t));
    var q = E.Times(pt, p % l);
    return E.Op(phi2, E.Op(t_phi, q));
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
        allBasis = EC.EllApFrobTrace(El.ArrModel.Select(e => e.Num).ToArray()).allBasis;

    GlobalStopWatch.AddLap();
    var N = EC.EllTateAlgorithm(EC.EllCoefs(curve)).N;
    var allList = EC.SmallPrimesList(N);
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
                Verif(Erl, phi, be1, be2, EC.Frob(be1), EC.Frob(be2), "phi");

                Console.WriteLine(new { phi2 });
                Verif(Erl, phi2, be1, be2, EC.Frob(be1, 2), EC.Frob(be2, 2), "phi^2");

                Console.WriteLine(new { add_phi2_pPy = add_phi2_p });
                Verif(Erl, add_phi2_p, be1, be2, Ep.Op(EC.Frob(be1, 2), Ep.Times(be1, p % l)),
                    Ep.Op(EC.Frob(be2, 2), Ep.Times(be2, p % l)), "phi^2 + [p]*P");
            }

            foreach (var t in l.SeqLazy())
            {
                var t_phi = Erl.FrobRl(Erl.Times(pt, -t));
                var eqFrob = Erl.Op(add_phi2_p, t_phi);
                if (log != LogLevel.Off)
                {
                    var (bg, be1, be2) = basis[l];
                    var Ep = Ell.ToGF(bg);
                    Verif(Erl, t_phi, be1, be2, EC.Frob(Ep.Times(be1, -t)), EC.Frob(Ep.Times(be2, -t)), "-[t]*phi");
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
            throw new();

        Console.WriteLine();
    }

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
    testLDivRing1();
    testLDivRing2();
}