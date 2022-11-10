using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class ConwayPolynoms
{
    static string ConwayWord(FpPolynom poly)
    {
        var p0 = poly.P;
        var digits = Enumerable.Repeat('0', $"{p0}".Length).Glue();
        var l = digits.Length;
        var fmt = $"{{0,{l}:{digits}}}";
        return poly.Coefs.Reverse().Select((e, i) => i % 2 == 0 ? e : IntExt.AmodP(p0 - e, p0)).Glue("", fmt);
    }

    static FpPolynom GetPoly(int p, int n)
    {
        var x = new FpPolynom(p);
        FpPolynom Poly_PpowD(int d) => x.Pow((int)Math.Pow(p, d)) - x;
        var ln = Enumerable.Range(0, p).MultiLoop(n).Select(l => l.Append(1).ToArray()).ToArray();

        var Px = Poly_PpowD(n);
        var lPolys = IntExt.Dividors(n).Select(Poly_PpowD).ToList();

        var irreductibles = new List<FpPolynom>();
        var zero = FpPolynom.Zero(p);

        foreach (var l in ln)
        {
            var tl = PolynomExt.TrimPoly(l);
            var lx = new FpPolynom(p, tl);
            if (lx.Equals(zero) || lx.Equals(x) || lx.CoefMax != 1)
                continue;

            if (Px.DivideBy(lx).rem.Equals(zero) &&
                lPolys.All(pi => !pi.DivideBy(lx).rem.Equals(zero)))
            {
                var gf0 = new GFp($"GF({lx})", ('x', p), tl);
                if (gf0.Count() == (int)Math.Pow(p, n) - 1)
                    irreductibles.Add(lx);
            }
        }

        var cnPoly = irreductibles.MinBy(ConwayWord);
        // foreach (var pi in irreductibles.OrderBy(ConwayWord))
        //     Console.WriteLine($"{pi} => {ConwayWord(pi)}");

        return cnPoly;
    }

    static void MyPoly(int p, int n)
    {
        var cnPoly = GetPoly(p, n);
        var g = PolynomExt.Get((int)Math.Pow(p, n));
        var cnPoly0 = new FpPolynom(p, g.coefs);

        var gf = new GFp($"GF({cnPoly})", ('x', p), cnPoly.Coefs);
        DisplayGroup.Head(Group.Create(gf));
        if (!cnPoly0.Equals(cnPoly))
        {
            Console.WriteLine($"My Poly     p={p} n={n} : {cnPoly} ### Conway Poly p={p} n={n} : {cnPoly0}");
            Console.WriteLine("############### Why ###############");
        }
        else
        {
            Console.WriteLine($"Conway Poly p={p} n={n} : {cnPoly}");
        }

        Console.WriteLine();
    }

    public static void AutomorphismFromPoly(int p, int n, bool verbose = true)
    {
        var cnPoly = GetPoly(p, n);
        Console.WriteLine($"Poly p={p} n={n} : {cnPoly}");
        var x = FpPolynom.Fpx(p);
        var q = (int)Math.Pow(p, n);
        var n1 = q - 1;
        var cn = new Cn(p);
        var gr = Product.GpGenerate(Enumerable.Repeat(cn, n).Cast<IGroup<ZnInt>>().ToArray());

        var list = Enumerable.Range(1, q).Select(i => x.Pow(i).DivideBy(cnPoly).rem.Coefs)
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

        Console.WriteLine("Check Distributivity {0}", distrib ? "Pass" : "Fail");
        Console.WriteLine("Check Fp-Space       {0}", fpspace ? "Pass" : "Fail");

        if (verbose)
            Console.WriteLine();
    }

    public static void Run()
    {
        MyPoly(2, 1);
        MyPoly(2, 2);
        MyPoly(2, 3);
        MyPoly(2, 4);
        MyPoly(2, 5);

        MyPoly(3, 1);
        MyPoly(3, 2);
        MyPoly(3, 2);
        MyPoly(3, 3);
        MyPoly(3, 4); // My Poly     p=3 n=4 : 2 + 2x + x^4 ### Conway Poly p=3 n=4 : 2 + 2x^3 + x^4
        MyPoly(3, 5);

        MyPoly(5, 1);
        MyPoly(5, 2);
        MyPoly(5, 3);
        MyPoly(5, 4); // My Poly     p=5 n=4 : 2 + 3x + x^2 + x^4 ### Conway Poly p=5 n=4 : 2 + 4x + 4x^2 + x^4

        MyPoly(7, 1);
        MyPoly(7, 2);
        MyPoly(7, 3); // My Poly     p=7 n=3 : 2 + 3x + x^3 ### Conway Poly p=7 n=3 : 4 + 6x^2 + x^3
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
        GlobalStopWatch.Time("F8", () => AutomorphismFromPoly(2, 3, verbose: false));
        GlobalStopWatch.Time("F9", () => AutomorphismFromPoly(3, 2, verbose: false));
        GlobalStopWatch.Time("F25", () => AutomorphismFromPoly(5, 2, verbose: false));
        GlobalStopWatch.Time("F27", () => AutomorphismFromPoly(3, 3, verbose: false));
        GlobalStopWatch.Time("F16", () => AutomorphismFromPoly(2, 4, verbose: false));
        GlobalStopWatch.Time("F125", () => AutomorphismFromPoly(5, 3, verbose: false));
        GlobalStopWatch.Time("F81", () => AutomorphismFromPoly(3, 4, verbose: false));
        GlobalStopWatch.Time("F32", () => AutomorphismFromPoly(2, 5, verbose: false));
    }
}