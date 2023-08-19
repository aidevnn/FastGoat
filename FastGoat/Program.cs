using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Matrix;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Polynoms;
using FastGoat.UserGroup.Words;
using FastGoat.UserGroup.Words.ToddCoxeter;
using static FastGoat.Commons.IntExt;
using static FastGoat.Commons.EnumerableExt;
using FastGoat.UserGroup.Padic;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
Ring.DisplayPolynomial = MonomDisplay.StarCaret;

HashSet<int>[] CyclotomicClasses(int q, int n)
{
    var cn = new Cn(n);
    var set = new HashSet<HashSet<int>>(new SetEquality<int>());
    var q0 = new ZnInt(n, q);
    foreach (var i0 in cn)
    {
        var Si = new HashSet<int>() { i0.K };
        var s = 1;
        while (true)
        {
            var qsi = q0.Pow(s) * i0;
            if (qsi.Equals(i0))
                break;

            Si.Add(qsi.K);
            ++s;
        }

        set.Add(Si.ToHashSet());
    }

    return set.ToArray();
}

(Fq fq, KPoly<ZnInt> g) GeneratorBCHmδ(int m, int delta)
{
    var fq = new Fq(2.Pow(m), 'a');
    var a = fq.X;
    var X = FG.KPoly('X', fq.X);
    // var nth = Group.MulGroup(fq.Name, a);
    // DisplayGroup.HeadElements(nth);
    var n = fq.Q - 1;
    var q = fq.P;

    var set = CyclotomicClasses(q, n);
    var set1 = set.Select(si => (si, si.Select(k => X - a.Pow(k)).Aggregate((s0, s1) => s0 * s1))).OrderBy(si => si.Item2).ToArray();
    set1.Println(si => $"[{si.si.Glue("; ")}] => {si.Item2}", $"{fq.FullName} nb cyclo classes = {set.Length}");

    var rg = (delta - 1).Range(1).Select(e0 => set1.First(e1 => e1.si.Contains(e0)).Item2).ToArray();
    var codeBCH = rg.Aggregate((e0, e1) => e0 * e1 / Ring.Gcd(e0, e1));
    rg.Println($"BCH(q:{q}, n:{n}, δ:{delta}) Code = {codeBCH}");
    Console.WriteLine();

    if (codeBCH.Coefs.Any(c => c.Poly.Degree != 0))
        throw new();

    return (fq, new('x', ZnInt.ZnZero(q), codeBCH.Coefs.Select(c => c[0]).ToArray()));
}

(Fq fq, KPoly<ZnInt> g) GeneratorBCHmd(int m, int d) => GeneratorBCHmδ(m, d + 1);

ZnInt[] RandMot(int k) => k.Range().Select(i => new ZnInt(2, Rng.Next(2))).ToArray();
KPoly<ZnInt> MotPoly(ZnInt[] mot) => new('x', ZnInt.ZnZero(2), mot);
KPoly<ZnInt> EncodeBCH(KPoly<ZnInt> m, KPoly<ZnInt> g) => m * g;

KPoly<EPoly<ZnInt>> Syndrom(KPoly<ZnInt> m, Fq fq, int delta)
{
    var b = fq.X;
    var X = FG.KPoly('X', b);
    return (delta - 1).Range(1).Select(j => m.Substitute(b.Pow(j)) * X.Pow(j)).Aggregate((e0, e1) => e0 + e1);
}

KPoly<ZnInt> Noise(KPoly<ZnInt> m, int delta)
{
    var d = m.Degree + 1;
    var t = (delta - 1) / 2;
    var m0 = new KPoly<ZnInt>(m.x, m.KZero, m.Coefs.ToArray());
    foreach (var j in t.Range().Select(i=>Rng.Next(d)).Distinct())
        m0.Coefs[j] = 1 - m0.Coefs[j];

    return m0;
}

(KPoly<EPoly<ZnInt>> u, KPoly<EPoly<ZnInt>> v) BerlekampMassey(KPoly<EPoly<ZnInt>> S, int delta)
{
    var t = (delta - 1) / 2;
    var (r0, u0, v0) = (S.One, S.Zero, S.X.Pow(2 * t));
    var (r1, u1, v1) = (S.Zero, S.One, S);
    while (!(v0.Degree >= t && v1.Degree < t))
    {
        var (quo, rem) = v0.Div(v1);
        // vi−1 = vi*qi + vi+1 puis les soustractions ri+1 = ri−1 − ri*qi et ui+1 = ui−1 − ui*qi
        (v0, v1) = (v1, rem);
        (r0, r1) = (r1, r0 - r1 * quo);
        (u0, u1) = (u1, u0 - u1 * quo);
    }

    return (u1, v1);
}

{
    // GAP examples 
    // CyclotomicClasses(2, 21).Println(si => $"[{si.Glue("; ")}]");
    // CyclotomicClasses(10, 21).Println(si => $"[{si.Glue("; ")}]");

    // Wikipedia examples
    // GeneratorBCHmd(4, 4);
    // GeneratorBCHmd(4, 6);
    // GeneratorBCHmd(4, 8);

    // Algebre Tome1, page 391-394
    var delta = 7;
    var (fq, g) = GeneratorBCHmδ(5, delta);
    var k = fq.Q - g.Degree;
    var t = (delta - 1) / 2;
    for (int i = 1; i <= 5; i++)
    {
        var m = EncodeBCH(MotPoly(RandMot(k)), g);
        var mn = Noise(m, delta);
        var S = Syndrom(mn, fq, delta);
        Console.WriteLine($"Mot{i}    : {m.Coefs.Glue().PadRight(fq.Q, '0')}");
        Console.WriteLine($"  Noise : {mn.Coefs.Glue().PadRight(fq.Q, '0')}");
        Console.WriteLine($"  Error : {(mn - m).Coefs.Glue().PadRight(fq.Q, '0')}");
        Console.WriteLine($"  Syndrom : {S}");
        var (u, v) = BerlekampMassey(S, delta);
        Console.WriteLine($"  BerlekampMassey : u = {u} and v = {v}");

        var x2t = S.X.Pow(2 * t);
        var v1 = (u * S).Div(x2t).rem;
        Console.WriteLine($"  Verif : {v.Equals(v1)}");
        if (!v.Equals(v1))
            throw new();
    }
}
