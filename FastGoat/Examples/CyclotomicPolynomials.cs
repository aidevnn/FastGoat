using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.Examples;

public static class CyclotomicPolynomials
{
    public static void First30CyclotomicsPolynomials()
    {
        for (int i = 1; i <= 30; i++)
        {
            Console.WriteLine($"Phi({i}) = {FG.CyclotomicPolynomial(i)}");
        }
    }

    public static void First20NthRootQ()
    {
        for (int n = 3; n <= 22; n++)
        {
            var unq = new NthRootQ(n);
            DisplayGroup.HeadElements(unq);
            DisplayGroup.AreIsomorphics(unq, new Cn(n)); // TODO fix fail for n = 2
            Console.WriteLine();
        }
    }

    public static void First20NthRootFq()
    {
        // for (int q = 9; q <= 9; q++)
        foreach (var q in new[] { 4, 8 })
        {
            if (IntExt.PrimesDecomposition(q).Distinct().Count() > 1)
                continue;

            for (int n = 45; n <= 45; n++)
            {
                if (IntExt.Gcd(q, n) != 1)
                    continue;

                var unq = new NthRootFq(n, q);
                DisplayGroup.HeadElements(unq);
                DisplayGroup.AreIsomorphics(unq, new Cn(n)); // TODO fix fail for n = 2
                Console.WriteLine();
            }
        }
    }
}