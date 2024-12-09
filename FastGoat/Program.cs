using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup;
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

// Learning with errors LWE encryption
// https://github.com/predrag3141/pslq_vs_lwe
bool RunTest(int nIn, int mIn, int maxQin, int O)
{
    foreach (var q in Primes10000.Where(q => q >= mIn && q <= maxQin))
    {
        var rgM = mIn.Range();
        Console.WriteLine($"############## q = {q,-3} n = {nIn,-3} m = {mIn,-3}");
        var arrA = Ring.Matrix(nIn, (mIn * nIn).Range().Select(_ => Rng.Next(q)).ToArray());
        var arrX = new int[mIn];
        var squaredNormTarget = (q * q) / 16.0;
        Console.WriteLine("A");
        Ring.DisplayMatrix(arrA);
        while (arrX.Sum(e => e * e) < squaredNormTarget)
        {
            var indexToModify = Rng.Next(mIn);
            arrX[indexToModify] += 2 * Rng.Next(2) - 1;
        }

        Console.WriteLine($"Private Key x:[{arrX.Glue(", ")}]");
        var rawU = nIn.Range().Select(i => mIn.Range().Sum(j => arrX[j] * arrA[i, j])).ToArray();

        var u = rawU.Select(e => AmodP(e, q)).ToArray();
        var rawUOverQ = u.Select((e, i) => (e - rawU[i]) / q).ToArray();
        Console.WriteLine($"rawU:     [{rawU.Glue(", ")}]");
        Console.WriteLine($"u:        [{u.Glue(", ")}]");
        Console.WriteLine($"rawUOverQ:[{rawUOverQ.Glue(", ")}]");

        var bs = (int)(0.5 + double.Pow(20, (mIn + nIn + 0.0) / nIn));
        Console.WriteLine($"base:{bs}");

        var v = (mIn + nIn + 1).Range().Select(_ => BigInteger.Zero).ToArray();
        var powerOfBase = BigInteger.One;
        foreach (var i in nIn.Range())
        {
            v[i] = q * powerOfBase;
            foreach (var j in mIn.Range())
                v[j + nIn] += powerOfBase * arrA[i, j];

            v[nIn + mIn] -= powerOfBase * u[i];
            powerOfBase *= bs;
        }

        Console.WriteLine("Input v to PSLQ");
        Console.WriteLine($"[{v.Glue(", ")}]");
        int[] w = [..rawUOverQ, ..arrX, 1];
        Console.WriteLine($"Expected output w of PSLQ:[{w.Glue(", ")}]");
        var vDotX = rgM.Select(i => v[nIn + i] * arrX[i]).Aggregate((a0, a1) => a0 + a1);
        Console.WriteLine($"<v,x> - v[{mIn}] = {vDotX - v[nIn + mIn]}");
        var vDotW = (mIn + nIn + 1).Range().Select(i => v[i] * w[i]).Aggregate((a0, a1) => a0 + a1);
        Console.WriteLine($"<v,w> = {vDotW}");
        if (vDotW != 0)
            throw new();

        var one = BigReal.BrOne(O);
        var gamma = 2 / BigReal.Sqrt(3 * one);
        BigInteger[] arrY;
        try
        {
            var vMat = v.Select(e => BigReal.FromBigInteger(e, O)).ToKMatrix();
            var arrY0 = PSLQM2<Dble, BigReal>.TwoLevelMultipair(vMat, gamma, vMat.N - 1);
            arrY = arrY0.Select(e => e.RoundEven.ToRational.Num).ToArray();
            if (arrY.Length == 0)
                continue;
        }
        catch (Exception e)
        {
            Console.WriteLine(e);
            continue;
        }

        Console.WriteLine("End PSLQM2 algorithm");
        Console.WriteLine("Possible Solution");
        Console.WriteLine($"[{arrY.Glue(", ")}]");

        if (arrY[mIn + nIn] == -1)
        {
            arrY = arrY.Select(e => -e).ToArray();
            Console.WriteLine("Get opposite");
            Console.WriteLine($"[{arrY.Glue(", ")}]");
        }

        Console.WriteLine();
        var sum = (mIn + nIn + 1).Range().Select(i => arrY[i] * v[i]).Aggregate((a0, a1) => a0 + a1);
        if (sum != 0)
        {
            Console.WriteLine("PSLQ found an incorrect solution");
            break;
        }

        var solutionIsCausal = true;
        foreach (var i in nIn.Range())
        {
            var yDotAi = mIn.Range().Select(j => arrY[nIn + j] * arrA[i, j]).Aggregate((a0, a1) => a0 + a1);
            if (BigInteger.Remainder(yDotAi - u[i], q) != 0)
            {
                solutionIsCausal = false;
                var ydotAi = (int)BigInteger.Remainder(yDotAi, q);
                Console.WriteLine($"<y,A[{i}]> = {AmodP(ydotAi, q)} != {u[i]} mod {q}");
            }
        }

        if (solutionIsCausal)
        {
            if (arrY[nIn + mIn] != 1)
            {
                Console.WriteLine(
                    $"PSLQ Found a solution [{arrY.Glue(", ")}] with coefficient for u {arrY[nIn + mIn]} != 1");
            }
            else
            {
                var diffSqNorm = mIn.Range().Select(i => BigInteger.Pow(arrY[i + nIn] - arrX[i], 2))
                    .Aggregate((a0, a1) => a0 + a1);
                if (diffSqNorm == 0)
                {
                    Console.WriteLine($"PSLQ found the private key !!!!");
                    Console.Beep();
                    return true;
                }
                else
                {
                    var norm2 = (double)mIn.Range().Select(i => BigInteger.Pow(arrY[i], 2))
                        .Aggregate((a0, a1) => a0 + a1);
                    var norm = double.Sqrt(norm2);
                    if (norm < q / 4.0)
                    {
                        Console.WriteLine($"PSLQ Found a causal solution with norm {norm:g4} y != x");
                    }
                    else
                    {
                        Console.WriteLine($"PSLQ Found a solution with norm {norm:g4} > q/4 = {q / 4.0}");
                    }
                }
            }
        }
        else
        {
            Console.WriteLine("Non causal");
        }
    }

    Console.WriteLine();

    return false;
}

{
    var O = 40; // digits for PSLQ
    var n = 3;
    var m = (int)(2 * n * double.Log2(n));
    var maxQ = 10 * m;

    for (int i = 0; i < 1000; i++)
    {
        if (RunTest(n, m, maxQ, O))
        {
            Console.WriteLine($"Attempt:{i + 1,4}/1000");
            break;
        }
    }
}