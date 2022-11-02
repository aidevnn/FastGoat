using System.Threading.Channels;
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
    static void MultiplicationTable<T>(ConcreteGroup<T> gr, bool verbose = true) where T : struct, IElt<T>
    {
        var primes = gr.GetGenerators().Select(x => IntExt.PrimesDecomposition(Group.Generate(gr, x).Count()))
            .ToArray();
        
        if (primes.All(x => x.Count() != 1) || primes.SelectMany(x => x).Distinct().Count() != 1)
            throw new GroupException(GroupExceptionType.GroupDef);
        
        var aut = Group.AutBase(gr);
        var n = gr.Count() - 1;
        var fb = Group.AllMorphisms(gr, gr, Group.MorphismType.Isomorphism).Select(aut.Create)
            .First(e => Group.Cycle(aut, e).Count == n);

        var zero = gr.Neutral();
        var one = gr.GetGenerators().Aggregate(gr.Op);
        var pows = Group.Cycle(aut, fb).ToDictionary(e => e.Key[one], e => e.Value);
    
        T Mul(T a1, T a2)
        {
            if (a1.Equals(gr.Neutral()) || a2.Equals(gr.Neutral()))
                return gr.Neutral();

            var i = pows[a1];
            var j = pows[a2];
            return aut.Times(fb, i + j)[one];
        }

        bool CheckDistributivity(T m, T a1, T a2)
        {
            var e0 = Mul(m, gr.Op(a1, a2));
            var e1 = gr.Op(Mul(m, a1), Mul(m, a2));
            return e0.Equals(e1);
        }

        Console.WriteLine("Multiplication Table of F{0}", gr.Count());
        var distrib = true;
        foreach (var e1 in gr)
        {
            foreach (var e2 in gr)
            {
                var e3 = Mul(e1, e2);
                if (distrib)
                    distrib &= gr.All(m => CheckDistributivity(m, e1, e2));

                if (verbose)
                    Console.WriteLine("{0} x {1} = {2}", e1, e2, e3);

            }

            if (verbose)
                Console.WriteLine();
        }

        Console.WriteLine("Check Distributivity {0}", distrib ? "Pass" : "Fail");

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
        MultiplicationTable(Product.Generate(c3, c3));
        MultiplicationTable(Product.Generate(c3, c3, c3)); // 70sec
        MultiplicationTable(Product.Generate(c5, c5)); // 32sec
        MultiplicationTable(Product.Generate(c2, c3)); // throw exception
    }

    public static void Bench()
    {
        var c2 = new Cn(2);
        var c3 = new Cn(3);
        var c5 = new Cn(5);
        
        GlobalStopWatch.Time("F8", () => MultiplicationTable(Product.Generate(c2, c2, c2), verbose: false));
        GlobalStopWatch.Time("F9", () => MultiplicationTable(Product.Generate(c3, c3), verbose: false));
        GlobalStopWatch.Time("F19", () => MultiplicationTable(new Cn(19), verbose: false));
        GlobalStopWatch.Time("F25", () => MultiplicationTable(Product.Generate(c5, c5), verbose: false));
        GlobalStopWatch.Time("F27", () => MultiplicationTable(Product.Generate(c3, c3, c3), verbose: false));
    }
}