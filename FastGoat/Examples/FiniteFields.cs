using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;
using FastGoat.UserGroup.Words;

namespace FastGoat.Examples;

public static class FiniteFields
{
    static Automorphism<T> Candidate<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        var primes = gr.GetGenerators().Select(x => IntExt.PrimesDecomposition(Group.Generate(gr, x).Count()))
            .ToArray();

        if (primes.All(x => x.Count() != 1) || primes.SelectMany(x => x).Distinct().Count() != 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        var aut = Group.AutBase(gr);
        var n = gr.Count() - 1;
        var fb = Group.AllMorphisms(gr, gr, Group.MorphismType.Isomorphism).Select(aut.Create)
            .First(e => Group.Cycle(aut, e).Count == n);

        return fb;
    }

    static void MultiplicationTable<T>(ConcreteGroup<T> gr, bool verbose = true) where T : struct, IElt<T>
    {
        var decomp = IntExt.PrimesDecomposition(gr.Count()).ToArray();
        var characteristic = decomp[0];
        var n = gr.Count() - 1;
        var zero = gr.Neutral();
        var one = gr.GetGenerators().Ascending().First(); // .Aggregate(gr.Op);

        var fb = Candidate(gr);
        var aut = fb.AutGroup;
        var elt2pow = Group.Cycle(aut, fb).ToDictionary(e => e.Key[one], e => e.Value);
        var pow2elt = elt2pow.ToDictionary(e => e.Value, e => e.Key);

        T Mul(T a1, T a2)
        {
            if (a1.Equals(zero) || a2.Equals(zero))
                return zero;

            var i = elt2pow[a1];
            var j = elt2pow[a2];

            // return aut.Times(fb, i + j)[one];
            var k = (i + j) % n;
            return pow2elt[k == 0 ? n : k];
        }

        bool CheckDistributivity(T m, T a1, T a2)
        {
            var e0 = Mul(m, gr.Op(a1, a2));
            var e1 = gr.Op(Mul(m, a1), Mul(m, a2));
            return e0.Equals(e1);
        }

        bool CheckFpSpace(int k, T a1, T a2)
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
        var fps = Enumerable.Range(0, characteristic).ToArray();
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

    static void MinPoly(int q, bool matrixForm = false)
    {
        var decomp = IntExt.PrimesDecomposition(q).ToArray();
        if (decomp.Distinct().Count() != 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        var dim = decomp.Length;
        var characteristic = decomp[0];

        var cn = new Cn(characteristic);
        var gr = Product.GpGenerate(Enumerable.Repeat(cn, dim).Cast<IGroup<ZnInt>>().ToArray());
        var n = gr.Count() - 1;

        var fb = Candidate(gr);
        var aut = fb.AutGroup;

        var gens = gr.GetGenerators().Ascending().ToArray();
        var one = gens[0];
        var zero = gr.Neutral();
        var elt2pow = Group.Cycle(aut, fb).ToDictionary(e => e.Key[one], e => e.Value);
        var pow2elt = elt2pow.ToDictionary(e => e.Value, e => e.Key);

        Ep<ZnInt> Mul(Ep<ZnInt> a1, Ep<ZnInt> a2)
        {
            if (a1.Equals(zero) || a2.Equals(zero))
                return zero;

            var i = elt2pow[a1];
            var j = elt2pow[a2];

            // return aut.Times(fb, i + j)[one];
            var k = (i + j) % n;
            return pow2elt[k == 0 ? n : k];
        }

        Ep<ZnInt> ExtMul(int k, Ep<ZnInt> e)
        {
            var ke = e.Ei.Select(e0 => cn.Times(e0, k)).ToArray();
            return new(ke);
        }

        List<Ep<ZnInt>> FqCycle(Ep<ZnInt> e)
        {
            List<Ep<ZnInt>> seq = new() { one };
            while (true)
            {
                var last = seq.Last();
                var e1 = Mul(e, last);
                seq.Add(e1);
                if (e1.Equals(one))
                    break;
            }

            return seq;
        }

        string Comb(Ep<ZnInt> e)
        {
            return e.Ei.Select((z, i) => (z.K, i)).Where(a => a.K != 0)
                .Select(a => a.K == 1 ? $"e[{a.i}]" : $"{a.K}*e[{a.i}]").Glue(" + ");
        }

        string PolyStr(int k, int m)
        {
            if (m == 0)
                return $"{k}";

            var ks = k == 1 ? "" : $"{k}";
            var xm = m == 1 ? "a" : $"a^{m}";
            return $"{ks}{xm}";
        }

        Console.WriteLine("Candidate for F{0}", gr.Count());
        var gFq = gr.Except(new[] { zero }).Ascending().First(e => FqCycle(e).Count() == n + 1);
        Console.WriteLine("Primitive element a   = {0}", gFq);
            
        var cycle = FqCycle(gFq).Take(dim + 1).ToArray();
        for (int j = 2; j < dim + 1; ++j)
            Console.WriteLine("                  a^{0} = {1}", j, cycle[j]);
            
        var vectors = Product.GpGenerate(Enumerable.Repeat(cn, dim + 1).Cast<IGroup<ZnInt>>().ToArray());
        var minPoly = vectors.Where(v =>
                !v.Equals(vectors.Neutral()) &&
                v.Ei.Reverse().Select((z, i) => ExtMul(z.K, cycle[i])).Aggregate(gr.Op).Equals(zero))
            .Ascending().First();

        Console.WriteLine("Minimal Polynomial  : {0} = 0",
            minPoly.Ei.Reverse().Select((z, i) => (z.K, i)).Where(a => a.K != 0).Select(a => PolyStr(a.K, a.i))
                .Reverse().Glue(" + "));

        if (matrixForm)
        {
            Console.WriteLine("   * * *");
            for (int i = 0; i < dim; ++i)
            {
                var e = gens[i];
                var fe = fb[e];
                Console.WriteLine("e[{0}] = {1} => {2} = f[{0}] = {3}", i, e, fe, Comb(fe));
            }
        }

        Console.WriteLine();
    }

    public static void Fp()
    {
        MultiplicationTable(new Cn(2));
        MultiplicationTable(new Cn(3));
        MultiplicationTable(new Cn(5));
        MultiplicationTable(new Cn(7));
        MultiplicationTable(new Cn(19));
        MultiplicationTable(new Cn(8)); // throw exception
    }

    public static void PolyCandidate()
    {
        var matrixForm = false;
        MinPoly(2, matrixForm);
        MinPoly(3, matrixForm);
        MinPoly(5, matrixForm);
        MinPoly(7, matrixForm);
        MinPoly(11, matrixForm);
        MinPoly(13, matrixForm);
        
        MinPoly(4, matrixForm);
        MinPoly(9, matrixForm);
        MinPoly(25, matrixForm);
        MinPoly(49, matrixForm);
        MinPoly(121, matrixForm);
        MinPoly(169, matrixForm);

        MinPoly(8, matrixForm);
        MinPoly(27, matrixForm);
        MinPoly(125, matrixForm);
        MinPoly(343, matrixForm);

        MinPoly(16, matrixForm);
        MinPoly(81, matrixForm);
    }

    public static void FpWhyNot()
    {
        var s3 = new Sn(3);
        var c3 = Group.Generate("C3", s3, s3[(1, 2, 3)]);
        MultiplicationTable(c3);

        var s5 = new Sn(5);
        var c5 = Group.Generate("C5", s5, s5[(1, 2, 3, 4, 5)]);
        MultiplicationTable(c5);

        var s11 = new Sn(11);
        var c11 = Group.Generate("C11", s11, s11.Cycle(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11));
        MultiplicationTable(c11);

        // Ambiguous notation
        MultiplicationTable(new WordGroup("a5"));
        MultiplicationTable(new WordGroup("a7"));
    }

    public static void Fq()
    {
        var c2 = new Cn(2);
        var c3 = new Cn(3);
        var c5 = new Cn(5);

        MultiplicationTable(Product.Generate(c2, c2));
        MultiplicationTable(Product.Generate(c2, c2, c2));
        MultiplicationTable(Product.Generate(c2, c2, c2, c2));
        MultiplicationTable(Product.Generate(c3, c3));
        MultiplicationTable(Product.Generate(c3, c3, c3));
        MultiplicationTable(Product.Generate(c5, c5));
        MultiplicationTable(Product.Generate(c2, c3)); // throw exception
    }

    public static void Bench()
    {
        var c2 = new Cn(2);
        var c3 = new Cn(3);
        var c5 = new Cn(5);

        MultiplicationTable(new Cn(19), verbose: false); // starting up

        GlobalStopWatch.Time("F19", () => MultiplicationTable(new Cn(19), verbose: false));
        GlobalStopWatch.Time("F8", () => MultiplicationTable(Product.Generate(c2, c2, c2), verbose: false));
        GlobalStopWatch.Time("F9", () => MultiplicationTable(Product.Generate(c3, c3), verbose: false));
        GlobalStopWatch.Time("F25", () => MultiplicationTable(Product.Generate(c5, c5), verbose: false));
        GlobalStopWatch.Time("F27", () => MultiplicationTable(Product.Generate(c3, c3, c3), verbose: false));
        GlobalStopWatch.Time("F16", () => MultiplicationTable(Product.Generate(c2, c2, c2, c2), verbose: false));
        GlobalStopWatch.Time("F125", () => MultiplicationTable(Product.Generate(c5, c5, c5), verbose: false));
        GlobalStopWatch.Time("F81", () => MultiplicationTable(Product.Generate(c3, c3, c3, c3), verbose: false));
        GlobalStopWatch.Time("F32", () => MultiplicationTable(Product.Generate(c2, c2, c2, c2, c2), verbose: false));
    }
}