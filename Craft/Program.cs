using System.Numerics;
using System.Text;
using Examples;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
using FastGoat.UserGroup.EllCurve;
using FastGoat.UserGroup.Floats;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;
using static FastGoat.Commons.IntExt;

//////////////////////////////////
//                              //
//    Work Table for crafting   //
//                              //
//////////////////////////////////

Console.WriteLine("Hello World");
RecomputeAllPrimesUpTo(5000000);

BigInteger RhoF(BigInteger x, BigInteger N) => BigInteger.Remainder(x * x + 1, N);

BigInteger RhoFactor(BigInteger N)
{
    BigInteger x1 = 2;
    var x2 = RhoF(x1, N);
    BigInteger d;
    var set = new HashSet<(BigInteger, BigInteger)>();
    set.Add((x1, x2));
    do
    {
        x1 = RhoF(x1, N);
        x2 = RhoF(RhoF(x2, N), N);
        d = GcdBigInt(AmodPbigint(x1 - x2, N), N);
        Console.WriteLine($"    d={d} x1={x1} x2={x2} N={N}");
        if (!set.Add((x1, x2)))
            break;
    } while (!(d > 1 && d < N));

    set.Clear();
    if (d == 1)
        d = N;
    
    return d;
}

IEnumerable<BigInteger> RhoPrime(BigInteger N)
{
    var N0 = N;
    while (N0 != 1)
    {
        var d = RhoFactor(N0);
        yield return d;
        N0 /= d;
    }
}

{
    for (BigInteger i = 3; i < 100; i++)
    {
        var facts = RhoPrime(i).Order().ToArray();
        Console.WriteLine($"i = {i} = {facts.Glue(" * ")} = {PrimesDecompositionBigInt(i).Glue(" * ")}");
        Console.WriteLine();
    }

    {
        var i = 144493;
        var facts = RhoPrime(i).Order().ToArray();
        Console.WriteLine($"i = {i} = {facts.Glue(" * ")} = {PrimesDecompositionBigInt(i).Glue(" * ")}");
        Console.WriteLine();
    }
    
    {
        BigInteger i = BigInteger.Pow(10, 6) + 3;
        var facts = RhoPrime(i).Order().ToArray();
        Console.WriteLine($"i = {i} = {facts.Glue(" * ")} = {PrimesDecompositionBigInt(i).Glue(" * ")}");
        Console.WriteLine();
    }
}