using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.LWE;
using FastGoat.UserGroup.Polynoms;
using Rq = FastGoat.Structures.VecSpace.KPoly<FastGoat.UserGroup.Integers.Rational>;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RngSeed(259128);
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

Rq InvRq(Rq a, Rq f, Rational t)
{
    var a0 = a.ResModSigned(f, t);
    var ia = Ring.FastBezout(a0, f).Item1.ResModSigned(f, t);
    var q = (a0 * ia).ResModSigned(f, t);
    if (q.Degree != 0 || q.IsZero())
        throw new($"a={a} not invertible");

    var iq = InvModPbezbigint(q[0].Num, t.Num);
    var iaf = (new Rational(iq) * ia).CoefsModSigned(t);
    return iaf;
}

(Rq[] invs, Rq[] facts) CrtBase(int N, Rational t)
{
    if (!int.IsPow2(N))
        throw new();
    if (!Primes10000.Contains((int)t.Num) || !t.Mod(N).IsOne())
        throw new();
    var n = N / 2;
    var t0 = (int)t.Num;
    var sols = SolveAll_k_pow_m_equal_one_mod_n_strict_long(t0, N)
        .Select(e => new Rational(e).Signed(t0))
        .Order()
        .ToArray();
    
    var pm = FG.QPoly().Pow(n) + 1;
    var facts = sols.Select(e => pm.X - e).ToArray();
    var invs = facts.Select(mi => (quo: pm / mi, mi)).Select(e => e.quo * InvRq(e.quo, e.mi, t)).ToArray();
    facts.Println($"pm = {pm}");
    Console.WriteLine();

    return (invs, facts);
}

Rq SolveCRTRq(Rq[] ais, Rq[] invs, Rq pm, Rational t)
{
    return ais.Zip(invs).Aggregate(pm.Zero, (acc, ej) => (acc + ej.First * ej.Second).ResModSigned(pm, t));
}

Rq[] InvCRTRq(Rq x, Rq[] mods, Rational t) => mods.Select(mod => x.ResModSigned(mod, t)).ToArray();

{
    RecomputeAllPrimesUpTo(1000000);
    var x = FG.QPoly();
    for (int i = 2; i < 9; ++i)
    {
        var N = 1 << i;
        var n = N / 2;
        var t0 = Primes10000.First(t0 => t0 % N == 1);
        var t = new Rational(t0);
        var pm = x.Pow(n) + 1;

        Console.WriteLine($"N = {N} pm = {pm} t = {t0}");
        var crtBase = CrtBase(N, t);
        
        for (int j = 0; j < 5; j++)
        {
            var m1 = RLWE.GenUnif(n, t);
            var crt1 = InvCRTRq(m1, crtBase.facts, t).ToArray();
            crt1.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m1 = {m1}");
            var rm1 = SolveCRTRq(crt1, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm1}");
            Console.WriteLine();
            
            var m2 = RLWE.GenUnif(n, t);
            var crt2 = InvCRTRq(m2, crtBase.facts, t).ToArray();
            crt2.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m2 = {m2}");
            var rm2 = SolveCRTRq(crt2, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm2}");
            Console.WriteLine();

            var m1m2 = (m1 * m2).ResModSigned(pm, t);
            var crt12a = InvCRTRq(m1m2, crtBase.facts, t).ToArray();
            crt12a.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m1m2a = {m1m2}");
            var rm1m2a = SolveCRTRq(crt12a, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm1m2a}");
            Console.WriteLine();

            var crt12b = crt1.Zip(crt2).Zip(crtBase.facts)
                .Select(e => (e.First.First * e.First.Second).ResModSigned(e.Second, t)).ToArray();
            crt12b.Zip(crtBase.facts).Println(e => $"m = {e.First,-6} mod ({e.Second}) mod {t}", $"m1m2b = {m1m2}");
            var rm1m2b = SolveCRTRq(crt12b, crtBase.invs, pm, t);
            Console.WriteLine($"  => {rm1m2b}");
            Console.WriteLine();
            
            if (!m1.Equals(rm1) || !m2.Equals(rm2) || !m1m2.Equals(rm1m2a) || !m1m2.Equals(rm1m2b))
                throw new("Fail");
        }
    }
}