using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;

namespace FastGoat.Examples;

public static class PolynomialFactorization
{
    private static Random rnd = new Random();

    static KPoly<ZnInt> ProdIrr(int p, int d)
    {
        var x = FG.ZPoly(p);
        return x.Pow(p.Pow(d)) - x;
    }

    static bool IsIrreductibleFp(KPoly<ZnInt> f)
    {
        if (f.Degree < 1)
            return true;

        var p = f.P;
        var n = f.Degree;
        var divs = IntExt.Dividors(n);
        return ProdIrr(p, n).Div(f).rem.IsZero() && divs.All(d => !ProdIrr(p, d).Div(f).rem.IsZero());
    }

    static KPoly<ZnInt> GpowP(KPoly<ZnInt> f)
    {
        var p = f.KZero.P;
        var test = f.Coefs.Select((c, i) => (c, i)).All(e => !e.c.IsZero() && e.i % p == 0);
        if (test)
        {
            var coefs = (f.Degree / p + 1).Range().Select(i => f.Coefs[p * i]).ToArray();
            return new KPoly<ZnInt>(f.x, f.KZero, coefs);
        }
        else
            return f;
    }

    static KPoly<ZnInt> RandPoly(int p, int n)
    {
        var coefs = (n + 1).Range().Select(i => new ZnInt(p, rnd.Next(p))).TrimSeq().ToArray();
        return new KPoly<ZnInt>('x', ZnInt.KZero(p), coefs);
    }

    public static void IrreductiblePolysF3()
    {
        var x = FG.ZPoly(3);

        var f0 = (x + 1) * (x + 2);
        Console.WriteLine($"{f0} {IsIrreductibleFp(f0)}; {f0.Pow(3)}");
        var f1 = x + 2;
        Console.WriteLine($"{f1} {IsIrreductibleFp(f1)}; {f1.Pow(3)}");
        var f2 = x.Pow(3) + 2 * x + 1;
        Console.WriteLine($"{f2} {IsIrreductibleFp(f2)}; {f2.Pow(3)}");
        var f3 = x.Pow(4) + 2 * x + 1;
        Console.WriteLine($"{f3} {IsIrreductibleFp(f3)}; {f3.Pow(3)}");
        var f4 = x.Pow(4) + 2 * x.Pow(3) + 2;
        Console.WriteLine($"{f4} {IsIrreductibleFp(f4)}; {f4.Pow(3)}");
        var f5 = (x.Pow(2) + 2 * x + 2) * (x + 1);
        Console.WriteLine($"{f5} {IsIrreductibleFp(f5)}; {f5.Pow(3)}");
    }

    public static void IrreductiblePolysF7()
    {
        var x = FG.ZPoly(7);

        var f0 = (x + 1) * (x + 2);
        Console.WriteLine($"{f0} {IsIrreductibleFp(f0)}");
        var f1 = x + 2;
        Console.WriteLine($"{f1} {IsIrreductibleFp(f1)}");
        var f2 = x.Pow(3) + 2 * x + 1;
        Console.WriteLine($"{f2} {IsIrreductibleFp(f2)}");
        var f3 = x.Pow(4) + 2 * x + 1;
        Console.WriteLine($"{f3} {IsIrreductibleFp(f3)}");
        var f4 = x.Pow(4) + 2 * x.Pow(3) + 2;
        Console.WriteLine($"{f4} {IsIrreductibleFp(f4)}");
        var f5 = (x.Pow(2) + 2 * x + 2) * (x + 1);
        Console.WriteLine($"{f5} {IsIrreductibleFp(f5)}");
    }

    public static void IrreductibleRandPolys()
    {
        var n = 5;
        var p = 3;
        for (int i = 0; i < 4 * n; i++)
        {
            var f = RandPoly(p, n);
            Console.WriteLine($"{{0,{-7 * (n + 1)}}} is irreductible : {{1}}", f, IsIrreductibleFp(f));
        }
    }

    static void Reduce(List<KPoly<ZnInt>> uks, KPoly<ZnInt> vk, List<KPoly<ZnInt>> hks)
    {
        Console.WriteLine(new { k = uks.Count, uk = uks.Last(), vk, hk = hks.Last() });
        var uk = uks.Last();
        var p = uk.KZero.P;
        if (uks.Count % p != 0)
        {
            var vk1 = Ring.Gcd(uk, vk).Monic;
            if (vk1.Degree == 0)
                return;

            uks.Add(uk / vk1);
            hks.Add(vk / vk1);
            Reduce(uks, vk1, hks);
        }
        else
        {
            uks.Add(uk / vk);
            Reduce(uks, vk, hks);
        }
    }

    static void SFmodP(KPoly<ZnInt> f)
    {
        var u = Ring.Gcd(f, f.Derivative).Monic;
        var v = f / u;
        var uks = new List<KPoly<ZnInt>>() { u };
        var hks = new List<KPoly<ZnInt>>() { f.One };
        Reduce(uks, v, hks);
    }

    public static void ReduceSFmodP()
    {
        var x = FG.ZPoly(3);
        var f = x.Pow(15) + 2 * x.Pow(14) + 2 * x.Pow(12) + x.Pow(11) + 2 * x.Pow(10) + 2 * x.Pow(8) + x.Pow(7) +
                2 * x.Pow(6) + 2 * x.Pow(4);

        SFmodP(f);
        var g = (x + 1).Pow(3) * (x.Pow(3) + 2 * x.Pow(2) + x + 2).Pow(2) * (x.Pow(2) + x + 2) * x.Pow(4);
        Console.WriteLine(f);
        Console.WriteLine(g);
    }

    static void SFQ(List<KPoly<Rational>> uk)
    {
        var f = uk.Last();
        var u = f / Ring.Gcd(f, f.Derivative).Monic;
        if (!f.Equals(u))
        {
            uk.Add(u);
            SFQ(uk);
        }
    }

    public static void ReduceSFQ()
    {
        var x = FG.QPoly();
        // f = x 12 + x 11 − x 9 − 2x 8 + x 5 + x 4
        var f = x.Pow(12) + x.Pow(11) - x.Pow(9) - 2 * x.Pow(8) + x.Pow(5) + x.Pow(4);
        Console.WriteLine(f);
        var uks = new List<KPoly<Rational>>() { f };
        SFQ(uks);
        Console.WriteLine("uks");
        Console.WriteLine(uks.Glue("\n"));
    }
}