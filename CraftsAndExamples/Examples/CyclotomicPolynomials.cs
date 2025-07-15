using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace CraftsAndExamples.Examples;

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
        for (int n = 2; n <= 22; n++)
        {
            var unq = new NthRootQ(n);
            DisplayGroup.HeadElements(unq);
            DisplayGroup.AreIsomorphics(unq, new Cn(n));
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
            DisplayGroup.AreIsomorphics(unq, new Cn(n));
            var primElt = unq.PrimitivesRoots().First();
            Console.WriteLine($"{unq} ~ F8({primElt}) ~ F2({primElt.KOne.X})({primElt})");
            Console.WriteLine($"    with {primElt.F} = 0 in {unq}");
            Console.WriteLine($"    and  {primElt.KOne.F} = 0 in F8");
            Console.WriteLine();
        }
    }

    public static void First20NthRootF9()
    {
        for (int n = 2; n <= 22; n++)
        {
            if (IntExt.Gcd(n, 3) != 1)
                continue;

            var unq = new NthRootFq(n, 9);
            DisplayGroup.HeadElements(unq);
            DisplayGroup.AreIsomorphics(unq, new Cn(n));
            var primElt = unq.PrimitivesRoots().First();
            Console.WriteLine($"{unq} ~ F9({primElt}) ~ F3({primElt.KOne.X})({primElt})");
            Console.WriteLine($"    with {primElt.F} = 0 in {unq}");
            Console.WriteLine($"    and  {primElt.KOne.F} = 0 in F9");
            Console.WriteLine();
        }
    }

    static void Fpm_Iso_FqZ(int p, int n, int d)
    {
        var m = n * d;
        var q1 = p.Pow(n);
        var q2 = p.Pow(m);
        var fq1 = Group.MulGroup($"F{q1}", FG.FqX(q1, 'α'));
        var fq2 = Group.MulGroup($"F{q2}", FG.FqX(q2, 'ζ'));
        var alpha = fq1.GetGenerators().First();
        var Z = FG.KPoly('ζ', alpha);
        var zeta = fq2.GetGenerators().First().F;
        var minPol = IntFactorisation.Firr(zeta.Substitute(Z), alpha).First();
        var fq3 = Group.MulGroup($"F{q1}({zeta.x})", FG.EPoly(minPol, 'y'));

        DisplayGroup.AreIsomorphics(fq2, fq3);
        Console.WriteLine($"   with {minPol} = 0 in {fq3}");
        Console.WriteLine();
    }

    // GF(p^nd) ~ GF(p^n)(ζ)
    public static void Examples_Fpm_Iso_FpnZ()
    {
        Fpm_Iso_FqZ(p: 2, n: 2, d: 3); // F64 ~ F4(ζ)
        Fpm_Iso_FqZ(p: 2, n: 3, d: 2); // F64 ~ F8(ζ)
        Fpm_Iso_FqZ(p: 2, n: 2, d: 4); // F256 ~ F4(ζ)
        Fpm_Iso_FqZ(p: 2, n: 4, d: 2); // F256 ~ F16(ζ)
        Fpm_Iso_FqZ(p: 3, n: 2, d: 3); // F729 ~ F9(ζ)
        Fpm_Iso_FqZ(p: 3, n: 3, d: 2); // F729 ~ F27(ζ)
    }
}