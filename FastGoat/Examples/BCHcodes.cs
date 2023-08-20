using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

// Algebre Tome1, page 387-397
public static class BCHcodes
{
    static HashSet<int>[] CyclotomicClasses(int q, int n)
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

    static (Fq fq, KPoly<ZnInt> g) GeneratorBCHmδ(int m, int delta)
    {
        var fq = new Fq(2.Pow(m), 'a');
        var a = fq.X;
        var X = FG.KPoly('X', fq.X);
        // var nth = Group.MulGroup(fq.Name, a);
        // DisplayGroup.HeadElements(nth);
        var n = fq.Q - 1;
        var q = fq.P;

        var set = CyclotomicClasses(q, n);
        var set1 = set.Select(si => (si, si.Select(k => X - a.Pow(k)).Aggregate((s0, s1) => s0 * s1))).OrderBy(si => si.Item2)
            .ToArray();
        set1.Println(si => $"[{si.si.Glue("; ")}] => {si.Item2}", $"{fq.FullName} nb cyclo classes = {set.Length}");

        var rg = (delta - 1).Range(1).Select(e0 => set1.First(e1 => e1.si.Contains(e0)).Item2).ToArray();
        var codeBCH = rg.Aggregate((e0, e1) => e0 * e1 / Ring.Gcd(e0, e1));
        rg.Println($"BCH(q:{q}, n:{n}, δ:{delta}) Code = {codeBCH}");
        Console.WriteLine();

        if (codeBCH.Coefs.Any(c => c.Poly.Degree != 0))
            throw new();

        return (fq, new('x', ZnInt.ZnZero(q), codeBCH.Coefs.Select(c => c[0]).ToArray()));
    }

    static (Fq fq, KPoly<ZnInt> g) GeneratorBCHmd(int m, int d) => GeneratorBCHmδ(m, d + 1);

    static ZnInt[] RandWord(int k) => k.Range().Select(i => new ZnInt(2, IntExt.Rng.Next(2))).ToArray();
    static KPoly<ZnInt> WordPoly(ZnInt[] mot) => new('x', ZnInt.ZnZero(2), mot);
    static KPoly<ZnInt> EncodeBCH(KPoly<ZnInt> m, KPoly<ZnInt> g) => m * g;

    static KPoly<EPoly<ZnInt>> Syndrom(KPoly<ZnInt> m, Fq fq, int delta)
    {
        var b = fq.X;
        var X = FG.KPoly('X', b);
        return (delta - 1).Range(1).Select(j => m.Substitute(b.Pow(j)) * X.Pow(j)).Aggregate((e0, e1) => e0 + e1);
    }

    static KPoly<ZnInt> Noise(KPoly<ZnInt> m, int delta)
    {
        var d = m.Degree + 1;
        var t = (delta - 1) / 2;
        var t0 = IntExt.Rng.Next(t);
        var m0 = new KPoly<ZnInt>(m.x, m.KZero, m.Coefs.ToArray());
        foreach (var j in t0.Range().Select(i => IntExt.Rng.Next(d)).Distinct())
            m0.Coefs[j] = 1 - m0.Coefs[j];

        return m0;
    }

    static (KPoly<EPoly<ZnInt>> u, KPoly<EPoly<ZnInt>> v) BerlekampMassey(KPoly<EPoly<ZnInt>> S, int delta)
    {
        var t = (delta - 1) / 2;
        var (r0, u0, v0) = (S.One, S.Zero, S.X.Pow(2 * t));
        var (r1, u1, v1) = (S.Zero, S.One, S);

        // vi−1 = vi*qi + vi+1 puis les soustractions ri+1 = ri−1 − ri*qi et ui+1 = ui−1 − ui*qi
        // jusqu’à obtenir deg vi < t et deg vi−1 >= t
        while (!(v1.Degree < t && v0.Degree >= t))
        {
            var (quo, rem) = v0.Div(v1);
            (v0, v1) = (v1, rem);
            (r0, r1) = (r1, r0 - r1 * quo);
            (u0, u1) = (u1, u0 - u1 * quo);
        }

        var c = u1.LT;
        return (u1 / c, v1 / c);
    }

    static KPoly<ZnInt> Error(KPoly<EPoly<ZnInt>> u, int n)
    {
        var b = u.KZero.X;
        var e = b.F.Zero.SubstituteChar('x');
        for (int i = 0; i < n - 1; i++)
        {
            var bi = b.Pow(i).Inv();
            if (u.Substitute(bi).IsZero())
                e += e.X.Pow(i);
        }

        return e;
    }

    static string Poly2Bin(KPoly<ZnInt> P, int n) => P.Coefs.Glue().PadRight(n, '0');

    static void CheckBCHmδ(int m, int delta, int nbTest)
    {
        var (fq, g) = GeneratorBCHmδ(m, delta);
        var n0 = fq.Q - g.Degree;
        Console.WriteLine($"BCH word {fq.Q}bit, max length {n0}");
        Console.WriteLine($"code {g}");
        Console.WriteLine();
        var t = (delta - 1) / 2;
        for (int k = 1; k <= nbTest; k++)
        {
            var send = EncodeBCH(WordPoly(RandWord(n0)), g);
            var received = Noise(send, delta);
            var err0 = received - send;
            var S = Syndrom(received, fq, delta);
            Console.WriteLine($"Syndrom : {S}");

            var b = S.KOne.X;
            var x = S.X;
            var I = err0.Coefs.Select((ei, i) => (ei, i)).Where(ei => !ei.ei.IsZero()).Select(ei => ei.i).ToArray();
            var u0 = I.Select(i => 1 - b.Pow(i) * x).Aggregate(x.One, (acc, xi) => acc * xi);
            var v0 = I.Select(i =>
                    b.Pow(i) * x * I.Except(new[] { i }).Select(j => 1 - b.Pow(j) * x).Aggregate(x.One, (acc, xj) => acc * xj))
                .Aggregate(x.Zero, (acc, xi) => acc + xi);

            var x2t = S.X.Pow(2 * t);
            var c = u0.LT;
            (u0, v0) = (u0 / c, v0 / c);
            Console.WriteLine($"  Direct          : u = {u0} and v = {v0}");
            var v2 = (u0 * S).Div(x2t).rem;
            if (!v0.Equals(v2))
                throw new();

            var (u, v) = BerlekampMassey(S, delta);
            Console.WriteLine($"  BerlekampMassey : u = {u} and v = {v}");

            var v1 = (u * S).Div(x2t).rem;
            Console.WriteLine($"  Verif Direct S : {v0.Equals(v2)}");
            Console.WriteLine($"  Verif B-M    S : {v.Equals(v1)}");
            if (!v.Equals(v1))
                throw new();
            Console.WriteLine($"  diff u : {u - u0}");
            Console.WriteLine($"  diff v : {v - v0}");

            var err1 = Error(u0, fq.Q);
            var err2 = Error(u, fq.Q);
            Console.WriteLine($"Send Word{k,-2}       : {Poly2Bin(send, fq.Q)}");
            Console.WriteLine($"  Receive         : {Poly2Bin(received, fq.Q)}");
            Console.WriteLine($"  Error Expected  : {Poly2Bin(err0, fq.Q)}");
            Console.WriteLine($"  Error Direct    : {Poly2Bin(err1, fq.Q)}");
            Console.WriteLine($"  Error Indirect  : {Poly2Bin(err2, fq.Q)}");

            Console.WriteLine();
        }
    }

    public static void BCH_Codes_Examples()
    {
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;

        CheckBCHmδ(m: 5, delta: 7, nbTest: 5); // 32bit
        CheckBCHmδ(m: 5, delta: 9, nbTest: 5); // 32bit
        
        CheckBCHmδ(m: 7, delta: 19, nbTest: 10); // 32bit
        CheckBCHmδ(m: 7, delta: 25, nbTest: 10); // 32bit
    }
    /*
        Syndrom : (a^6 + a^5 + a^4 + a^2)*X^24 + (a^3 + a + 1)*X^23 + (a^6 + a^5 + a^4 + a^3 + a^2 + a)*X^22 + (a^6 + a^5 + a^3 + a + 1)*X^21 + (a^5 + a^4)*X^20 + (a^4 + a^2 + 1)*X^19 + (a^6 + a^4 + a^3 + a^2 + 1)*X^18 + (a^6 + a^2 + 1)*X^17 + (a^5 + a^4)*X^16 + (a^5 + a^4 + a^3 + a)*X^15 + (a^6 + a^5 + a^4 + a^3 + a^2)*X^14 + (a^5 + a^3 + a)*X^13 + (a^6 + a^2 + a)*X^12 + (a^6 + a^5 + a^4)*X^11 + (a^6 + a^3 + a^2)*X^10 + (a^5 + a^3 + a + 1)*X^9 + (a^6 + a^3 + a^2)*X^8 + (a^6 + a^5 + a)*X^7 + (a^4 + a^3)*X^6 + (a^5 + a^3 + a^2 + a)*X^5 + (a^5 + a^3 + a^2 + a)*X^4 + a^5*X^3 + (a^6 + a^5 + a^4 + a^3 + a^2)*X^2 + (a^6 + a^5 + a)*X
          Direct          : u = X^9 + (a^6 + a^5 + a^3 + a^2 + a + 1)*X^8 + (a^4 + a)*X^7 + a^2*X^6 + (a^6 + a^5 + a^4 + 1)*X^5 + (a^6 + a^5 + a^2 + a + 1)*X^4 + (a^6 + a^2 + a)*X^3 + a^5*X^2 + X + a^5 + a^4 + a^3 + a and v = X^9 + (a^4 + a)*X^7 + (a^6 + a^5 + a^4 + 1)*X^5 + (a^6 + a^2 + a)*X^3 + X
          BerlekampMassey : u = X^9 + (a^6 + a^5 + a^3 + a^2 + a + 1)*X^8 + (a^4 + a)*X^7 + a^2*X^6 + (a^6 + a^5 + a^4 + 1)*X^5 + (a^6 + a^5 + a^2 + a + 1)*X^4 + (a^6 + a^2 + a)*X^3 + a^5*X^2 + X + a^5 + a^4 + a^3 + a and v = X^9 + (a^4 + a)*X^7 + (a^6 + a^5 + a^4 + 1)*X^5 + (a^6 + a^2 + a)*X^3 + X
          Verif Direct S : True
          Verif B-M    S : True
          diff u : 0
          diff v : 0
        Send Word5        : 01000100101111011111101010001011001100100110110100111010100000001111100101011100110110110000101011110000010001101101111010111010
          Receive         : 01000100111111011111101010011111001100100011100100111010100000001111100101011100110010110100101011110001010001101101111010111010
          Error Expected  : 00000000010000000000000000010100000000000101010000000000000000000000000000000000000100000100000000000001000000000000000000000000
          Error Direct    : 00000000010000000000000000010100000000000101010000000000000000000000000000000000000100000100000000000001000000000000000000000000
          Error Indirect  : 00000000010000000000000000010100000000000101010000000000000000000000000000000000000100000100000000000001000000000000000000000000
     */
}