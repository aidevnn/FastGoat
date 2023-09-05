using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class ConwayPolynoms
{
    static Dictionary<int, Dictionary<int, KPoly<ZnInt>>> allCnPolys = new();
    
    static KPoly<ZnInt> GetPoly(int p, int n)
    {
        var x = FG.ZPoly(p);
        var xPow = new Dictionary<int, KPoly<ZnInt>>();
        var xp = x;
        for (int k = 1; k <= n; k++)
        {
            var xp0 = xp.One;
            for (int j = 0; j < p; j++)
                xp0 *= xp;

            xp = xp0;
            xPow[k] = xp - x;
        }

        KPoly<ZnInt> Poly_PpowD(int d) => xPow[d];// x.Pow(p.Pow(d)) - x;
        var seq = EnumerableExt.MultiLoop(n.Range().OrderDescending()
            .Select(i => p.Range().Select(j => j == 0 ? x.Zero : x.Pow(i) * ((n - i) % 2 == 0 ? j : p - j))));

        var Px = Poly_PpowD(n);
        var lPolys = IntExt.Dividors(n).Select(m => (m, l: Poly_PpowD(m))).ToList();

        var cnPoly = x.Zero;
        foreach (var lt in seq)
        {
            var lx = x.Pow(n) + lt.Aggregate(x.Zero, (sum, xi) => xi + sum);
            if (lx.IsZero() || lx.Equals(x) || lx[0].IsZero() || !Px.Div(lx).rem.IsZero())
                continue;

            var a = FG.EPoly(lx);
            var pn = p.Pow(n) - 1;
            var aPow = new Dictionary<int, EPoly<ZnInt>>() { [0] = a.One };
            var acc = a.One;
            var i = 0;
            do
            {
                ++i;
                acc *= a;
                aPow[i] = acc;
            } while (!acc.Equals(a.One));

            if (i != pn)
                continue;

            if (lPolys.All(pi =>
                    !pi.l.Div(lx).rem.IsZero() &&
                    (n == 1 || allCnPolys[p][pi.m].Substitute(aPow[pn / (p.Pow(pi.m) - 1)]).IsZero())))
            {
                cnPoly = lx;
                break;
            }
        }

        if (!allCnPolys.ContainsKey(p))
            allCnPolys[p] = new() { [1] = cnPoly };
        else
            allCnPolys[p][n] = cnPoly;

        return cnPoly;
    }

    static void MyPoly(int p, int n)
    {
        var x = FG.ZPoly(p);
        var pn = p.Pow(n);
        var g = PolynomExt.GetConwayPoly(pn);
        var cnPoly0 = g.coefs.Select((k, i) => k * x.Pow(i)).Aggregate(x.Zero, (acc, xi) => acc + xi);
        var cnPoly1 = GetPoly(p, n);
        Console.WriteLine($"Conway Poly {$"p={p} n={n} |GF({cnPoly0})| = {pn - 1}", -80} Exact:{cnPoly1.Equals(cnPoly0)}");
        // var gf = new GFp($"GF({cnPoly})", FG.EPoly(FG.KPoly(ZnInt.KZero(p), 'x', cnPoly.Coefs)));
        // DisplayGroup.Head(Group.Generate(gf));
    }

    public static void AutomorphismFromPoly(int p, int n, bool verbose = true)
    {
        var cnPoly = GetPoly(p, n);
        if (verbose)
            Console.WriteLine($"Poly p={p} n={n} : {cnPoly}");
        
        var x = FG.ZPoly(p);
        var q = (int)Math.Pow(p, n);
        var n1 = q - 1;
        var cn = new Cn(p);
        var gr = Product.GpGenerate(Enumerable.Repeat(cn, n).Cast<IGroup<ZnInt>>().ToArray());

        var list = Enumerable.Range(1, q).Select(i => x.Pow(i).Div(cnPoly).rem.Coefs.Select(k => k.K).ToArray())
            .Select(c => c.Length < n ? c.Concat(new int[n - c.Length]).ToArray() : c).ToArray();
        var listEp = list.Select(c => Product.Ep(c.Select(i => cn[i]).ToArray())).ToArray();
        var pmap = listEp.SkipLast(1).Select((e, i) => (e, listEp[(i + 1) % q])).ToDictionary(e => e.e, e => e.Item2);

        var fbMap = Group.AutomorphismMap(gr, pmap);
        var aut = Group.AutBase(gr);
        var fb = aut.Create(fbMap);

        var fbCycle = Group.Cycle(aut, fb);
        var zero = gr.Neutral();
        var one = gr.GetGenerators().Ascending().First();
        var elt2pow = fbCycle.ToDictionary(e => e.Key[one], e => e.Value);
        var pow2elt = elt2pow.ToDictionary(e => e.Value, e => e.Key);

        Ep<ZnInt> Mul(Ep<ZnInt> a1, Ep<ZnInt> a2)
        {
            if (a1.Equals(zero) || a2.Equals(zero))
                return zero;

            var i = elt2pow[a1];
            var j = elt2pow[a2];

            // return aut.Times(fb, i + j)[one];
            var k = (i + j) % n1;
            return pow2elt[k == 0 ? n1 : k];
        }

        bool CheckDistributivity(Ep<ZnInt> m, Ep<ZnInt> a1, Ep<ZnInt> a2)
        {
            var e0 = Mul(m, gr.Op(a1, a2));
            var e1 = gr.Op(Mul(m, a1), Mul(m, a2));
            return e0.Equals(e1);
        }

        bool CheckFpSpace(int k, Ep<ZnInt> a1, Ep<ZnInt> a2)
        {
            var ka1 = gr.Times(a1, k);
            var ka2 = gr.Times(a2, k);
            var e0 = Mul(a1, a2);
            var ke0 = gr.Times(e0, k);
            var e1 = Mul(a1, ka2);
            var e2 = Mul(ka1, a2);
            return ke0.Equals(e1) && ke0.Equals(e2);
        }

        if (verbose)
            Console.WriteLine("Multiplication Table of F{0}", gr.Count());
        
        var distrib = true;
        var fpspace = true;
        var fps = Enumerable.Range(0, p).ToArray();
        foreach (var e1 in gr)
        {
            foreach (var e2 in gr)
            {
                var e3 = Mul(e1, e2);
                if (distrib)
                    distrib &= gr.All(m => CheckDistributivity(m, e1, e2));

                if (fpspace)
                    fpspace &= fps.All(k => CheckFpSpace(k, e1, e2));

                if (verbose)
                    Console.WriteLine("{0} x {1} = {2}", e1, e2, e3);
            }

            if (verbose)
                Console.WriteLine();
        }

        if (verbose)
        {
            Console.WriteLine("Check Distributivity {0}", distrib ? "Pass" : "Fail");
            Console.WriteLine("Check Fp-Space       {0}", fpspace ? "Pass" : "Fail");
            Console.WriteLine();
        }
    }

    // Generate Fq polynomial for q < 1024, and check validity 
    public static void RunTestFq()
    {
        var nb = 1024;
        foreach (var p in IntExt.Primes10000.Where(p => p * p < nb))
        {
            var mx = (int)(Double.Log(nb) / Double.Log(p));
            for (int n = 1; n <= mx; n++)
                MyPoly(p, n);
        }
    }

    public static void FastAutomorphism()
    {
        AutomorphismFromPoly(2, 1);
        AutomorphismFromPoly(3, 1);
        AutomorphismFromPoly(5, 1);
        AutomorphismFromPoly(7, 1);
        
        AutomorphismFromPoly(2, 2);
        AutomorphismFromPoly(3, 2);

        AutomorphismFromPoly(2, 3);
        AutomorphismFromPoly(3, 3);
    }

    public static void Bench()
    {
        var dico = new Dictionary<int, int>() { [2] = 7, [3] = 4, [5] = 3 };
        foreach (var (p, n) in dico)
        {
            foreach (var m in n.Range(1))
                GlobalStopWatch.Bench(4, $"F{p.Pow(m)}", () => AutomorphismFromPoly(p, m, verbose: false));
        }
    }
}