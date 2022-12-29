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

    public static void First20NthRootF8()
    {
        for (int n = 3; n <= 22; n++)
        {
            if (IntExt.Gcd(n, 2) != 1)
                continue;

            var unq = new NthRootFq(n, 8);
            DisplayGroup.HeadElements(unq);
            DisplayGroup.AreIsomorphics(unq, new Cn(n)); // TODO fix fail for n = 2
            Console.WriteLine();
        }
    }
}