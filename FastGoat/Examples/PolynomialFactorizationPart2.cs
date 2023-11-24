using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Padic;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.UserGroup.Polynoms.IntFactorisation;

namespace FastGoat.Examples;

public static class PolynomialFactorizationPart2
{
    public static void IrreductibleFactorizationZ()
    {
        Console.WriteLine();
        Console.WriteLine("Irreductible Factorization in Z[X]");
        Console.WriteLine();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var X = FG.QPoly('X');

        // AECF example 21.2, page 387
        FirrZ(X.Pow(4) - 1, details: true);
        FirrZ(X.Pow(15) - 1, details: true);
        FirrZ(X.Pow(11) - 1, details: true);
        FirrZ(X.Pow(8) - 40 * X.Pow(6) + 352 * X.Pow(4) - 960 * X.Pow(2) + 576, details: true);
        FirrZ(X.Pow(12) - 50 * X.Pow(10) + 753 * X.Pow(8) - 4520 * X.Pow(6) + 10528 * X.Pow(4) - 6720 * X.Pow(2) + 576, details: true);

        FirrZ((X + 3) * (X - 5) * (X + 11) * (X - 17), details: true);
        FirrZ(X.Pow(6) + 2 * X.Pow(4) - 1, details: true);
        FirrZ((X.Pow(3) + 3 * X.Pow(2) + -2) * (X.Pow(3) + -3 * X.Pow(2) + 2), details: true);
        FirrZ((X.Pow(3) + 3 * X.Pow(2) + -2) * (X.Pow(3) + -5 * X.Pow(2) + 2 * X + -4), details: true);

        FirrZ(X.Pow(6) + 5 * X.Pow(5) + 3 * X.Pow(4) + -7 * X.Pow(3) + -3 * X.Pow(2) + 7 * X + -2, details: true);

        FirrZ(X.Pow(12) + 5 * X.Pow(11) + -202 * X.Pow(10) + -155 * X.Pow(9) + 11626 * X.Pow(8) + -37275 * X.Pow(7)
              + -33479 * X.Pow(6) + 547100 * X.Pow(5) + -560012 * X.Pow(4) + -616520 * X.Pow(3) + 351876 * X.Pow(2)
              + 146520 * X + -56160, details: true);

        Console.WriteLine("Random Z[X] Polynomials");
        Console.WriteLine();
        for (int j = 0; j < 20; j++)
        {
            var amp = IntExt.Rng.Next(2, 9);
            var n = 12 + IntExt.Rng.Next(13);
            var degrees = IntExt.Partitions32[n].Where(l => l.All(i => i != 1) && l.Count > 1)
                .OrderBy(i => IntExt.Rng.NextDouble())
                .FirstOrDefault(new[] { 2, 3, 4 }.ToList())
                .ToArray();

            var polys = degrees.Select(ni => PolynomialFactorization.RandPolySep(Rational.KZero(), amp, ni, monic: false)).ToArray();
            var f0 = polys.Aggregate((a, b) => a * b);
            var f = new KPoly<Rational>('X', f0.KZero, f0.Coefs);
            if (Ring.Discriminant(f).IsZero()) // Separable polynomial
                continue;

            Console.WriteLine($"Random factors : [{polys.Glue(" ;  ")}]");
            FirrZ(f, details: true);

            /***
             *  Examples of outputs
             *
                f = X^6 + 12*X^5 + 38*X^4 + 8*X^3 + -9*X^2 + 73*X + -42
                Disc(f) = -95855189967891 = -1^1 * 3^20 * 37^1 * 743^1
                Prime P = 11; Sigma = 5; P^o=161051 Nu=12360.95
                f = X^6 + X^5 + 5*X^4 + -3*X^3 + 2*X^2 + -4*X + 2 mod 11
                Fact(f) = (X + 2)*(X + -4)*(X + -2)*(X^3 + 5*X^2 + -4*X + -4) mod 11
                Fact(f) = (X + 2)*(X^2 + 5*X + -3)*(X^3 + 5*X^2 + -4*X + 7) in Z[X]

                f = X^10 + 9*X^9 + 27*X^8 + 55*X^7 + -29*X^6 + -123*X^5 + -291*X^4 + 370*X^3 + 144*X^2 + -57*X + -126
                Disc(f) = 4406013110764543812064828211149450560 = 2^6 * 3^6 * 5^1 * 19^1 * 23^2 * 47^2 * 53^1 * 271^2
                #### Prime 7 and Sigma 8 wont work ####
                #### Prime 7 and Sigma 9 wont work ####
                Prime P = 13; Sigma = 6; P^o=4826809 Nu=1256602.80
                f = X^10 + -4*X^9 + X^8 + 3*X^7 + -3*X^6 + -6*X^5 + -5*X^4 + 6*X^3 + X^2 + -5*X + 4 mod 13
                Fact(f) = (X + 2)*(X + 3)*(X^2 + 3*X + 6)*(X^2 + -2*X + -4)*(X^4 + 3*X^3 + 2*X^2 + -5*X + -4) mod 13
                Fact(f) = (X^2 + 3*X + 6)*(X^2 + 5*X + -7)*(X^6 + X^5 + 5*X^4 + -8*X^3 + -2*X^2 + 2*X + 3) in Z[X]

             */
        }
    }

    public static void IrreductibleFactorizationLLL()
    {
        Console.WriteLine();
        Console.WriteLine("Irreductible Factorization in Z[X]");
        Console.WriteLine();
        Ring.DisplayPolynomial = MonomDisplay.StarCaret;
        var X = FG.QPoly('X');

        VanHoeijFactorization(X.Pow(4) - 1);
        VanHoeijFactorization(X.Pow(15) - 1);
        VanHoeijFactorization(X.Pow(11) - 1);
        VanHoeijFactorization(X.Pow(8) - 40 * X.Pow(6) + 352 * X.Pow(4) - 960 * X.Pow(2) + 576);
        VanHoeijFactorization(X.Pow(12) - 50 * X.Pow(10) + 753 * X.Pow(8) - 4520 * X.Pow(6) + 10528 * X.Pow(4) - 6720 * X.Pow(2) +
                              576);

        VanHoeijFactorization((X + 3) * (X - 5) * (X + 11) * (X - 17));
        VanHoeijFactorization(X.Pow(6) + 2 * X.Pow(4) - 1);
        VanHoeijFactorization((X.Pow(3) + 3 * X.Pow(2) + -2) * (X.Pow(3) + -3 * X.Pow(2) + 2));
        VanHoeijFactorization((X.Pow(3) + 3 * X.Pow(2) + -2) * (X.Pow(3) + -5 * X.Pow(2) + 2 * X + -4));

        VanHoeijFactorization(X.Pow(6) + 5 * X.Pow(5) + 3 * X.Pow(4) + -7 * X.Pow(3) + -3 * X.Pow(2) + 7 * X + -2);

        VanHoeijFactorization(X.Pow(12) + 5 * X.Pow(11) + -202 * X.Pow(10) + -155 * X.Pow(9) + 11626 * X.Pow(8) + -37275 * X.Pow(7)
                              + -33479 * X.Pow(6) + 547100 * X.Pow(5) + -560012 * X.Pow(4) + -616520 * X.Pow(3) + 351876 * X.Pow(2)
                              + 146520 * X + -56160);

        Console.WriteLine("Random Z[X] Polynomials");
        Console.WriteLine();
        for (int j = 0; j < 20; j++)
        {
            var amp = IntExt.Rng.Next(2, 29);
            var n = 2 + IntExt.Rng.Next(20);
            var degrees = IntExt.Partitions32[n].Where(l => l.All(i => i != 1) && l.Count > 1)
                .OrderBy(i => IntExt.Rng.NextDouble())
                .FirstOrDefault(new int[] { 2, 3, 4 }.ToList())
                .ToArray();

            var polys = degrees.Select(ni => PolynomialFactorization.RandPolySep(Rational.KZero(), amp, ni, monic: false)).ToArray();
            var f0 = polys.Aggregate((a, b) => a * b);
            var f = new KPoly<Rational>('X', f0.KZero, f0.Coefs);
            if (Ring.Discriminant(f).IsZero()) // Separable polynomial
                continue;

            Console.WriteLine($"Random factors : [{polys.Glue(" ;  ")}]");

            VanHoeijFactorization(f);

            /***
             *  Examples of outputs
             *
                f = X^6 + -10*X^5 + 28*X^4 + -20*X^3 + -65*X^2 + 50*X + 16
                Disc(f) = -1553297679360000 ~ -1^1 * 2^13 * 3^7 * 5^4 * 7^2 * 19^1 * 149^1
                Prime P = 11; Tau = 5; Theta = 661328432.0240374 Sigma = 10; Nu = 5.291502622129181
                f = X^6 + X^5 +  6*X^4 +  2*X^3 + X^2 +  6*X +  5 mod 11
                Fact(f) = (X^2 +  3*X +  6)*(X^2 +  6*X +  3)*(X +  4)*(X + 10) mod 11
                LLL Combinaisons
                    [0, 1, 0, 0]
                    [1, 0, 1, 0]
                    [0, 0, 0, 1]
                Fact(f) = (X^2 + -5*X + -8)*(X^3 + -4*X^2 + 7*X + 2)*(X + -1) in Z[X]
                f = Fact(f) : True

                f = X^9 + -16*X^8 + 68*X^7 + -38*X^6 + 68*X^5 + -404*X^4 + 1297*X^3 + -1372*X^2 + -1704*X + 1560
                Disc(f) = 171242376893905022721623419223750303809536 ~ 2^22 * 3^6 * 7^6 * 11^1 * 19^1 * 23^2 * 29^1 * 37^2 * 43^2 * 67^2 * 151^2 * 757^2
                Prime P = 5; Tau = 11; Theta = 8490538852469.95 Sigma = 22; Nu = 10.828203913853857
                f = X^9 + 4*X^8 + 3*X^7 + 2*X^6 + 3*X^5 + X^4 + 2*X^3 + 3*X^2 + X mod 5
                Fact(f) = (X + 2)*(X^2 + 3*X + 4)*(X + 4)*(X^2 + X + 1)*(X + 3)*(X + 1)*(X) mod 5
                LLL Combinaisons
                    [0, 1, 0, 0, 0, 0, 0]
                    [1, 0, 1, 0, 0, 0, 0]
                    [0, 0, 0, 1, 0, 0, 0]
                    [0, 0, 0, 0, 0, 1, 0]
                    [0, 0, 0, 0, 1, 0, 1]
                Fact(f) = (X^2 + -2*X + 4)*(X^2 + -9*X + 13)*(X^2 + -9*X + 6)*(X + 1)*(X^2 + 3*X + 5) in Z[X]
                f = Fact(f) : True

                f = X^11 + -17*X^10 + -16*X^9 + 767*X^8 + -1560*X^7 + 1397*X^6 + -491*X^5 + -2910*X^4 + 2799*X^3 + -1911*X^2 + -1685*X + -182
                Disc(f) = 72570018092906449167601721141392286400505133081053425861717760 ~ 2^8 * 5^1 * 7^1 * 11^3 * 23^2 * 113^2 * 139^2 * 173^2 * 233^1 * 2833^2 * 3449^1
                SearchVanHoeij P = X^11 + -17*X^10 + -16*X^9 + 767*X^8 + -1560*X^7 + 1397*X^6 + -491*X^5 + -2910*X^4 + 2799*X^3 + -1911*X^2 + -1685*X + -182 Tau = 17 Theta = 3402784769104.1543 BoundSigma = 2.4561005912838476E+175 MaxSigma = 368 => Nb = 5
                #### Prime 3 Sigma 34 = 2*tau wont work ####
                #### Prime 3 Sigma 51 = 3*tau wont work ####
                #### Prime 3 Sigma 68 = 4*tau wont work ####
                #### Prime 3 Sigma 85 = 5*tau wont work ####
                #### Prime 3 wont work ####
                SearchVanHoeij P = X^11 + -17*X^10 + -16*X^9 + 767*X^8 + -1560*X^7 + 1397*X^6 + -491*X^5 + -2910*X^4 + 2799*X^3 + -1911*X^2 + -1685*X + -182 Tau = 8 Theta = 212188261668464.22 BoundSigma = 1.3614502344021254E+195 MaxSigma = 176 => Nb = 5
                Prime P = 13; Tau = 8; Theta = 212188261668464.22 Sigma = 16; Nu = 10.246950765959598
                f = X^11 +  9*X^10 + 10*X^9 +  6*X^6 +  3*X^5 +  2*X^4 +  4*X^3 +  5*X mod 13
                Fact(f) = (X^6 +  9*X^5 +  5*X^4 +  7*X^3 + 10*X^2 + 11)*(X +  5)*(X + 10)*(X +  8)*(X +  3)*(X) mod 13
                LLL Combinaisons
                    [1, 0, 1, 0, 0, 0]
                    [0, 0, 0, 1, 1, 0]
                    [0, 1, 0, 0, 0, 1]
                Fact(f) = (X^7 + -7*X^6 + 4*X^5 + -8*X^4 + -11*X^3 + 9*X^2 + -15*X + -7)*(X^2 + -15*X + -2)*(X^2 + 5*X + -13) in Z[X]
                f = Fact(f) : True

             */
        }
    }
}