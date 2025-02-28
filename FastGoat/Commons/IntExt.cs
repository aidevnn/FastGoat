using System.Numerics;

namespace FastGoat.Commons;

/// <summary>
/// This static class provides extension methods for the <see cref="int"/> type.
/// </summary>
public static class IntExt
{
    /// <summary>
    /// Initializes the <see cref="IntExt"/> class.
    /// </summary>
    static IntExt()
    {
        Primes10000 = new(AllPrimes(100000));
        MaxPrime = Primes10000.Max();
        Partitions32 = IntPartitions(32);

        var comp = Comparer<int[]>.Create((a, b) => a.SequenceCompareTo(b));
        AllPermutations = new Dictionary<int, int[][]>
        {
            [0] = new[] { Array.Empty<int>() }
        };
        for (int k = 1; k <= NbPermutations; ++k)
        {
            int k0 = k;
            AllPermutations[k] = Enumerable.Range(0, k)
                .SelectMany(c =>
                    AllPermutations[k - 1].Select(l => l.Take(c).Append(k0).Concat(l.Skip(c)).ToArray()))
                .OrderBy(a => a, comp).ToArray();
        }

        SolveSquareInt = SquareSums();
        Rng = new();
    }

    /// <summary>
    /// Gets or sets the Random Rng.
    /// </summary>
    public static Random Rng { get; private set; }

    /// <summary>
    /// Gets a Random sign.
    /// </summary>
    public static int RngSign => (2 * Rng.Next(2) - 1);

    /// <summary>
    /// This method sets the seed value for the random number generator.
    /// </summary>
    /// <param name="seed">The seed value to be used for the random number generator.</param>
    public static void RngSeed(int seed)
    {
        Rng = new(seed);
    }

    public static Dictionary<int, Dictionary<int, int[][]>> SolveSquareInt { get; }

    public static long MaxPrime { get; private set; }
    /// <summary>
    /// All Primes less than 10000
    /// </summary>
    public static List<int> Primes10000 { get; private set; }

    /// <summary>
    ///  Partitions of an integers until 32. 
    /// </summary>
    public static Dictionary<int, List<List<int>>> Partitions32;


    /// <summary>
    /// Returns a sequence of prime numbers up to a given number.
    /// </summary>
    /// <param name="n">The maximum number to be included in the sequence.</param>
    /// <returns>A sequence of prime numbers up to the given number.</returns>
    static IEnumerable<int> AllPrimes(int n)
    {
        var primes = new Queue<int>();
        for (int i = 2; i <= n; ++i)
        {
            var isprime = true;
            foreach (var p in primes)
            {
                if (p * p > i) break;
                if (i % p == 0)
                {
                    isprime = false;
                    break;
                }
            }

            if (isprime) primes.Enqueue(i);
        }

        return primes;
    }

    /// <summary>
    /// Recomputing AllPrimes less than a given bound
    /// </summary>
    /// <param name="n">The maximum number to be included in the sequence.</param>
    public static void RecomputeAllPrimesUpTo(int n)
    {
        Console.WriteLine($"Recompute All Primes Up To {n}");
        Primes10000 = new(AllPrimes(n));
        MaxPrime = Primes10000.Max();
        Console.WriteLine("done");
        Console.WriteLine();
    }

    /// <summary>
    /// Primality test of a given integer. 
    /// </summary>
    /// <param name="n">The integer to test.</param>
    /// <returns>A boolean, true if the given integer is prime or false otherwise.</returns>
    public static bool IsPrime(BigInteger n)
    {
        if (n <= 1)
            return false;
        if (n < MaxPrime)
            return Primes10000.Contains((int)n);
        if (n > MaxPrime * MaxPrime)
            throw new($"n={n} is bigger than maxPrime^2={MaxPrime * MaxPrime}");

        return n > 1 && Primes10000.Select(e => new BigInteger(e)).Where(e => e * e <= n).All(e => n % e != 0);
    }

    /// <summary>
    /// Calculates the prime decomposition of a given Integer. 
    /// </summary>
    /// <param name="n">The Integer to decompose.</param>
    /// <returns>An IEnumerable of ints containing the prime factors of the given Integer.</returns>
    public static IEnumerable<int> PrimesDecomposition(int n)
    {
        var n0 = n;
        foreach (var p in Primes10000.Where(p => p <= n))
        {
            while (n0 % p == 0)
            {
                yield return p;
                n0 /= p;
            }

            if (n0 == 1) break;
        }

        if (n0 != 1)
        {
            var ns = double.Sqrt(n0);
            if (MaxPrime > ns)
                yield return n0;
            else
                throw new($"n0={n0} is bigger than maxPrime^2={MaxPrime * MaxPrime}");
        }
    }

    /// <summary>
    /// Calculates the prime decomposition of a given BigInteger. 
    /// </summary>
    /// <param name="n">The BigInteger to decompose.</param>
    /// <returns>An IEnumerable of ints containing the prime factors of the given BigInteger.</returns>
    public static IEnumerable<int> PrimesDecompositionBigInt(BigInteger n)
    {
        var n0 = n;
        foreach (var p in Primes10000.Where(p => p <= n))
        {
            while (n0 % p == 0)
            {
                yield return p;
                n0 /= p;
            }

            if (n0 == 1) yield break;
        }

        if (n0 != 1)
        {
            var ns = double.Sqrt((double)n0);
            if (MaxPrime > ns && n0 < int.MaxValue)
                yield return (int)n0;
            else if (MaxPrime < ns)
                throw new($"n0={n0} is bigger than maxPrime^2={MaxPrime * MaxPrime}");
            else
                throw new($"n0={n0} is prime but bigger than intMax={int.MaxValue}");
        }
    }

    /// <summary>
    /// Calculates the multiplicity of a given factor in a given number and the remaining factor, a=p^m * r 
    /// </summary>
    /// <param name="a">The number to decompose.</param>
    /// <param name="p">The factor.</param>
    /// <returns>A tuple of integer, the multiplicity of factor p and the remaining factor.</returns>
    public static (int mul, BigInteger rem) FactorMultiplicity(int p, BigInteger a)
    {
        var mul = 0;
        var rem = a;
        while (rem % p == 0)
        {
            ++mul;
            rem /= p;
        }

        return (mul, rem);
    }

    /// <summary>
    /// This method returns a dictionary containing the prime factors of a given number.
    /// </summary>
    /// <param name="n">The number to be factored.</param>
    /// <returns>A Dictionary containing the prime factors of the given number.</returns>
    public static Dictionary<int, int> PrimesDec(BigInteger n)
    {
        var dec = PrimesDecompositionBigInt(BigInteger.Abs(n));
        var dico = dec.GroupBy(e => e).ToDictionary(e => e.Key, e => e.Count());
        if (n < 0)
            dico[-1] = 1;
        return dico;
    }

    /// <summary>
    /// Partitions of an integers until N
    /// Return a dictionary where each key is mapped to a collection
    /// </summary>
    /// <param name="n">int</param>
    /// <returns>A dictionary of all partitions</returns>
    /// <example>5 => {5}, {4,1}, {3,2}, {3,1,1}, {2,2,1}, {2,1,1,1}, {1,1,1,1,1}</example>
    public static Dictionary<int, List<List<int>>> IntPartitions(int n)
    {
        List<List<int>> all = new();
        Queue<List<int>> q = new();
        q.Enqueue(new());

        while (q.Count != 0)
        {
            var l0 = q.Dequeue();
            var s = l0.Count == 0 ? 1 : l0.First();
            for (int k = s; k <= n; ++k)
            {
                var l1 = l0.ToList();
                l1.Insert(0, k);
                var sumL1 = l1.Sum();
                if (sumL1 <= n)
                {
                    q.Enqueue(l1);
                    all.Add(l1);
                }
                else
                    break;
            }
        }

        return all.GroupBy(l0 => l0.Sum()).ToDictionary(a => a.Key, b => b.ToList());
    }

    /// <summary>
    /// Map of solutions of equation
    /// <list type="bullet">
    /// <item><term>a^2=n</term></item>
    /// <item><term>a^2+b^2=n</term></item>
    /// <item><term>a^2+b^2+c^2=n</term></item>
    /// </list>
    /// </summary>
    /// <param name="maxSum"></param>
    /// <param name="maxList"></param>
    /// <returns>A dictionary of solutions</returns>
    public static Dictionary<int, Dictionary<int, int[][]>> SquareSums(int maxSum = 128, int maxList = 10)
    {
        var maxN = (int)Double.Sqrt(maxSum);
        var allSums = (maxN - 1).Range(2).ToDictionary(i => i * i, i => new List<List<int>>() { new() { i } });
        while (true)
        {
            var tmp = allSums.ToDictionary(a => a.Key, a => a.Value.Select(b => b.ToList()).ToList());
            foreach (var (sum, list) in allSums)
            {
                foreach (var l1 in list.Where(l => l.Count < maxList))
                {
                    var s = l1.Last();
                    for (int s1 = s; s1 <= maxN; s1++)
                    {
                        var sum2 = s1 * s1 + sum;
                        if (sum2 <= maxSum)
                        {
                            var l2 = l1.Append(s1).ToList();
                            if (!tmp.ContainsKey(sum2))
                                tmp[sum2] = new List<List<int>>() { l2 };
                            else if (tmp[sum2].All(l => !l.SequenceEqual(l2)))
                                tmp[sum2].Add(l2);
                        }
                    }
                }
            }

            if (allSums.Sum(k => k.Value.Count) != tmp.Sum(k => k.Value.Count))
                allSums = tmp;
            else
                break;
        }

        var squareSums = allSums.SelectMany(e => e.Value.Select(l => (sum: e.Key, l.Count, l)))
            .OrderBy(e => e.Count)
            .GroupBy(e => e.Count)
            .ToDictionary(
                e => e.Key,
                e => e.OrderBy(f => f.sum).GroupBy(f => f.sum).ToDictionary(
                    f => f.Key,
                    f => f.Select(g => g.l.ToArray())
                        .OrderBy(l => l, Comparer<int[]>.Create((la, lb) => la.SequenceCompareTo(lb)))
                        .ToArray()));

        return squareSums;
    }

    /// <summary>
    /// Represents the maximum of precomputed permutations.
    /// </summary>
    public const int NbPermutations = 8;

    private static Dictionary<int, int[][]> AllPermutations { get; }

    /// <summary>
    /// Insertion and Iteration Based Permutation Generation
    /// Generates all possible permutations of a set of numbers. 
    /// </summary>
    /// <param name="n">The number of elements in the set.</param>
    /// <returns>An enumerable of enumerables containing all possible permutations.</returns>
    public static IEnumerable<IEnumerable<int>> YieldAllPermutations(int n)
    {
        if (n == 0)
            yield return Enumerable.Empty<int>();
        else
            foreach (var perm in YieldAllPermutations(n - 1))
                for (int i = 0; i < n; i++)
                    yield return perm.InsertAt(i, n);
    }

    /// <summary>
    /// Generates all combinations of k elements from a set of n elements.
    /// </summary>
    /// <param name="k">The number of elements in each combination.</param>
    /// <param name="n">The total number of elements in the set.</param>
    /// <returns>An enumerable containing all possible combinations.</returns>
    public static IEnumerable<IEnumerable<bool>> YieldCombsKinN(int k, int n)
    {
        if (k == 0)
            yield return Enumerable.Repeat(false, n);
        else if (k == n)
            yield return Enumerable.Repeat(true, n);
        else
        {
            foreach (var list in YieldCombsKinN(k, n - 1))
                yield return list.InsertAt(n - 1, false);

            foreach (var list in YieldCombsKinN(k - 1, n - 1))
                yield return list.InsertAt(n - 1, true);
        }
    }

    /// <summary>
    /// Generates all possible combinations of boolean values for a given number of elements.
    /// </summary>
    /// <param name="n">The number of elements.</param>
    /// <returns>An enumerable of enumerables containing all possible combinations of boolean values.</returns>
    public static IEnumerable<IEnumerable<bool>> YieldAllCombs(int n)
    {
        for (int k = 0; k <= n; k++)
        {
            foreach (var list in YieldCombsKinN(k, n))
            {
                yield return list;
            }
        }
    }

    /// <summary>
    /// Generates all possible combinations of boolean values for a given number of elements from m to n.
    /// </summary>
    /// <param name="m">The minimal number of true elements.</param>
    /// <param name="n">The number of elements.</param>
    /// <returns>An enumerable of enumerables containing all possible combinations of boolean values.</returns>
    public static IEnumerable<IEnumerable<bool>> YieldAllCombsFromMtoN(int m, int n)
    {
        for (int k = m; k <= n; k++)
        {
            foreach (var list in YieldCombsKinN(k, n))
            {
                yield return list;
            }
        }
    }

    /// <summary>
    /// Generates all possible combinations of boolean values for a given number of elements from 1 to m.
    /// </summary>
    /// <param name="m">The maximal number of true elements.</param>
    /// <param name="n">The number of elements.</param>
    /// <returns>An enumerable of enumerables containing all possible combinations of boolean values.</returns>
    public static IEnumerable<IEnumerable<bool>> YieldAllCombsToMofN(int m, int n)
    {
        for (int k = 0; k <= m; k++)
        {
            foreach (var list in YieldCombsKinN(k, n))
            {
                yield return list;
            }
        }
    }

    /// <summary>
    /// Generates all possible combinations of boolean values for a given number of elements.
    /// </summary>
    /// <param name="n">The number of elements.</param>
    /// <returns>An enumerable of enumerables containing all possible combinations of boolean values.</returns>
    public static IEnumerable<IEnumerable<bool>> YieldAllCombsBinary(int n)
    {
        IEnumerable<bool> IntToBinary(int n0, int p)
        {
            var p0 = p;
            for (int i = 0; i < n0; i++)
            {
                yield return (p0 & 1) == 1;
                p0 >>= 1;
            }
        }

        var mx = 1 << n;
        for (int i = 0; i < mx; i++)
        {
            yield return IntToBinary(n, i);
        }
    }

    /// <summary>
    /// Generates all possible combinations of boolean values for a given number of elements.
    /// </summary>
    /// <param name="n">The number of elements.</param>
    /// <returns>An enumerable of enumerables containing all possible combinations of boolean values.</returns>
    public static IEnumerable<IEnumerable<bool>> YieldAllCombinations(int n)
    {
        return new[] { false, true }.MultiLoop(n);
    }

    /// <summary>
    /// Calculates the value of n to the power of m. 
    /// </summary>
    /// <param name="n">The base number.</param>
    /// <param name="m">The exponent.</param>
    /// <returns>The result of n to the power of m.</returns>
    public static int Pow(this int n, int m) => (int)Math.Pow(n, m);

    /// <summary>
    /// Calculates the greatest common divisor of two integers. 
    /// </summary>
    /// <param name="a">The first integer.</param>
    /// <param name="b">The second integer.</param>
    /// <returns>The greatest common divisor of the two integers.</returns>
    public static int Gcd(int a, int b)
    {
        if (a == 0 && b == 0)
            throw new($"Gcd(0,0) dont exists");

        if (a * b == 0)
            return (a + b) < 0 ? -(a + b) : a + b;

        var q = a / b;
        var r = a - b * q;
        return Gcd(b, r);
    }

    /// <summary>
    /// Calculates the greatest common divisor of two integers. 
    /// </summary>
    /// <param name="a">The first integer.</param>
    /// <param name="b">The second integer.</param>
    /// <returns>The greatest common divisor of the two integers.</returns>
    public static long GcdLong(long a, long b)
    {
        if (a == 0 && b == 0)
            throw new($"Gcd(0,0) dont exists");

        if (a * b == 0)
            return (a + b) < 0 ? -(a + b) : a + b;

        var q = a / b;
        var r = a - b * q;
        return GcdLong(b, r);
    }

    /// <summary>
    /// Calculates the greatest common divisor of a given array of integers.
    /// </summary>
    /// <param name="arr">The array of integers to calculate the GCD for.</param>
    /// <returns>The greatest common divisor of the given array.</returns>
    public static int Gcd(int[] arr)
    {
        if (!arr.Any())
            return 0;

        return Gcd(arr.First(), Gcd(arr.Skip(1).ToArray()));
    }

    /// <summary>
    /// Calculates the greatest common divisor (GCD) of two BigIntegers.
    /// </summary>
    /// <param name="a">The first BigInteger.</param>
    /// <param name="b">The second BigInteger.</param>
    /// <returns>The GCD of the two BigIntegers.</returns>
    public static BigInteger GcdBigInt(BigInteger a, BigInteger b)
    {
        if (b == 0)
            return a < 0 ? -a : a;

        return GcdBigInt(b, BigInteger.Remainder(a, b));
    }

    /// <summary>
    /// Calculates the greatest common divisor (GCD) of an array of BigIntegers.
    /// </summary>
    /// <param name="arr">An array of BigIntegers.</param>
    /// <returns>The GCD of the given array.</returns>
    public static BigInteger GcdBigInt(BigInteger[] arr)
    {
        if (!arr.Any())
            return 0;

        return GcdBigInt(arr.First(), GcdBigInt(arr.Skip(1).ToArray()));
    }

    /// <summary>
    /// Calculates the least common multiple of two BigIntegers.
    /// </summary>
    /// <param name="a">The first BigInteger.</param>
    /// <param name="b">The second BigInteger.</param>
    /// <returns>The least common multiple of the two BigIntegers.</returns>
    public static BigInteger LcmBigInt(BigInteger a, BigInteger b) => (a * b) / BigInteger.GreatestCommonDivisor(a, b);

    /// <summary>
    /// Calculates the least common multiple of an array of BigIntegers.
    /// </summary>
    /// <param name="arr">The array of BigIntegers.</param>
    /// <returns>The least common multiple of the array.</returns>
    public static BigInteger LcmBigInt(BigInteger[] arr)
    {
        if (!arr.Any())
            return 1;

        return LcmBigInt(arr.First(), LcmBigInt(arr.Skip(1).ToArray()));
    }

    // wikipedia
    /// <summary>
    /// Computes the Bezout coefficients and greatest common divisor of two BigIntegers.
    /// </summary>
    /// <param name="a">The first BigInteger.</param>
    /// <param name="b">The second BigInteger.</param>
    /// <returns>A tuple containing the Bezout coefficients (Xa, Xb) and the greatest common divisor (Gcd).</returns>
    public static (BigInteger Xa, BigInteger Xb, BigInteger Gcd) BezoutBigInt(BigInteger a, BigInteger b)
    {
        var (r0, x0, y0, r1, x1, y1) = (a, BigInteger.One, BigInteger.Zero, b, BigInteger.Zero, BigInteger.One);
        while (!r1.IsZero)
        {
            var q = r0 / r1;
            (r0, x0, y0, r1, x1, y1) = (r1, x1, y1, r0 - q * r1, x0 - q * x1, y0 - q * y1);
        }

        if (r0 > 0)
            return (x0, y0, r0);

        return (-x0, -y0, -r0);
    }

    /// <summary>
    /// Calculates the Bezout coefficients for two integers a and b. 
    /// </summary>
    /// <param name="a">The first integer.</param>
    /// <param name="b">The second integer.</param>
    /// <returns>A tuple containing the Bezout coefficients for a and b.</returns>
    public static (long x, long y) BezoutLong(long a, long b)
    {
        if (b == 0)
        {
            var x = a < 0 ? -1 : 1;
            return (x, 0);
        }

        var q = a / b;
        var (x0, y0) = BezoutLong(b, a - b * q);
        return (y0, x0 - y0 * q);
    }

    /// <summary>
    /// Calculates the Bezout coefficients for two integers a and b. 
    /// </summary>
    /// <param name="a">The first integer.</param>
    /// <param name="b">The second integer.</param>
    /// <returns>A tuple containing the Bezout coefficients for a and b.</returns>
    public static (int x, int y) Bezout(int a, int b)
    {
        if (b == 0)
        {
            var x = a < 0 ? -1 : 1;
            return (x, 0);
        }

        var q = a / b;
        var (x0, y0) = Bezout(b, a - b * q);
        return (y0, x0 - y0 * q);
    }

    /// <summary>
    /// Calculates the Bezout coefficients (x, y) for two integers a and b. 
    /// </summary>
    /// <param name="a">The first integer.</param>
    /// <param name="b">The second integer.</param>
    /// <returns>A tuple containing the Bezout coefficients (x, y).</returns>
    public static (int x, int y) BezoutVerbose(int a, int b)
    {
        if (b == 0)
        {
            var x = a < 0 ? -1 : 1;
            Console.WriteLine($"End  {new { a, b, x, y = 0 }}");
            return (x, 0);
        }

        // gcd = a.x1 + b.y1
        // gcd = b.x0 + r.y0
        // gcd = b.x0 + (a-bq).y0
        // gcd = a.y0 + b(x0-q.y0)

        var q = a / b;
        var r = a - b * q;
        var (x0, y0) = BezoutVerbose(b, r);
        Console.WriteLine($"Step {new { a, b, q, r, x = y0, y = x0 - y0 * q }}");
        return (y0, x0 - y0 * q);
    }

    /// <summary>
    /// Calculates the result of raising a number to a power and then modulo it with a given number. 
    /// </summary>
    /// <param name="a">The base number.</param>
    /// <param name="exp">The exponent.</param>
    /// <param name="mod">The modulus.</param>
    /// <returns>The result of the calculation.</returns>
    public static int PowMod(int a, int exp, int mod)
    {
        if (exp < 0)
        {
            var ai = InvModPbez(a, mod);
            return PowMod(ai, -exp, mod);
        }
        
        var (r, a0, e0) = (1, AmodP(a, mod), exp);
        while (e0 > 0)
        {
            if (e0 % 2 == 1)
                r = (r * a0) % mod;
            e0 >>= 1;
            a0 = (a0 * a0) % mod;
        }

        return AmodP(r, mod);
    }

    /// <summary>
    /// Calculates the result of raising a number to a power and then modulo it with a given number. 
    /// </summary>
    /// <param name="a">The base number.</param>
    /// <param name="exp">The exponent.</param>
    /// <param name="mod">The modulus.</param>
    /// <returns>The result of the calculation.</returns>
    public static long PowModLong(long a, long exp, long mod)
    {
        if (exp < 0)
        {
            var ai = InvModPbezlong(a, mod);
            return PowModLong(ai, -exp, mod);
        }

        var (r, a0, e0) = ((long)1, a % mod, exp);
        while (e0 > 0)
        {
            if (e0 % 2 == 1)
                r = (r * a0) % mod;
            e0 >>= 1;
            a0 = (a0 * a0) % mod;
        }

        return AmodPlong(r, mod);
    }
    
    /// <summary>
    /// Calculates the result of raising a number to a power and then modulo it with a given number. 
    /// </summary>
    /// <param name="a">The base number.</param>
    /// <param name="exp">The exponent.</param>
    /// <param name="mod">The modulus.</param>
    /// <returns>The result of the calculation.</returns>
    public static BigInteger PowModBigint(BigInteger a, BigInteger exp, BigInteger mod)
    {
        if (exp < 0)
        {
            var ai = InvModPbezbigint(a, mod);
            return PowModBigint(ai, -exp, mod);
        }

        return AmodPbigint(BigInteger.ModPow(a, exp, mod), mod);
    }

    /// <summary>
    /// Calculates the Legendre-Jacobi of two given numbers (m|n). 
    /// </summary>
    /// <param name="m">The numerator.</param>
    /// <param name="n">The denominator.</param>
    /// <returns>The result of the calculation +/- 1.</returns>
    public static int LegendreJacobi(int m, int n) => PowMod(m, (n - 1) / 2, n);

    /// <summary>
    /// Calculates the Legendre-Jacobi of two given numbers (m|n). 
    /// </summary>
    /// <param name="m">The numerator.</param>
    /// <param name="n">The denominator.</param>
    /// <returns>The result of the calculation +/- 1.</returns>
    public static long LegendreJacobiLong(long m, long n) => PowModLong(m, (n - 1) / 2, n);

    /// <summary>
    /// Calculates the Legendre-Jacobi of two given numbers (m|n). 
    /// </summary>
    /// <param name="m">The numerator.</param>
    /// <param name="n">The denominator.</param>
    /// <returns>The result of the calculation +/- 1.</returns>
    public static BigInteger LegendreJacobiBigint(BigInteger m, BigInteger n) => PowModBigint(m, (n - 1) / 2, n);

    /// <summary>
    /// Solves the system of congruences using the Chinese Remainder Theorem (CRT).
    /// </summary>
    /// <param name="a">Array of remainders.</param>
    /// <param name="m">Array of moduli.</param>
    /// <returns>A tuple containing the solution x and the product of all moduli mod.</returns>
    /// <exception cref="ArgumentException">Thrown when the arrays have different lengths or the moduli are not pairwise coprime.</exception>
    public static (int x, int mod) CRTint(int[] a, int[] m)
    {
        if (m.Length < 0 || m.Length != a.Length)
            throw new ArgumentException("Array dimension");

        if (m.Any(mi => mi < 2) || m.Grid2D().Where(e => e.t1 != e.t2).Any(e => Gcd(e.t1, e.t2) != 1))
            throw new ArgumentException("Modulus");

        var A = a.Zip(m).Select(e => AmodP(e.First, e.Second)).ToArray();
        var mod = m.Aggregate((mi, mj) => mi * mj);
        var x = 0;
        foreach (var (mi, ai) in m.Zip(A))
        {
            var quo = mod / mi;
            var inv = InvModPbez(quo, mi);
            x = (x + ai * quo * inv) % mod;
        }

        return (x % mod, mod);
    }

    /// <summary>
    /// Solves the system of congruences using the Chinese Remainder Theorem (CRT).
    /// </summary>
    /// <param name="a">Array of remainders.</param>
    /// <param name="m">Array of moduli.</param>
    /// <returns>A tuple containing the solution x and the product of all moduli mod.</returns>
    /// <exception cref="ArgumentException">Thrown when the arrays have different lengths or the moduli are not pairwise coprime.</exception>
    public static (long x, long mod) CRTlong(long[] a, long[] m)
    {
        if (m.Length < 0 || m.Length != a.Length)
            throw new("array dimension");

        if (m.Any(mi => mi < 2) || m.Grid2D().Where(e => e.t1 != e.t2).Any(e => GcdLong(e.t1, e.t2) != 1))
            throw new("modulus");

        var A = a.Zip(m).Select(e => AmodPlong(e.First, e.Second)).ToArray();
        var mod = m.Aggregate((mi, mj) => mi * mj);
        long x = 0;
        foreach (var (mi, ai) in m.Zip(A))
        {
            var quo = mod / mi;
            var inv = InvModPbezlong(quo, mi);
            x = (x + ai * quo * inv) % mod;
        }

        return (x, mod);
    }

    /// <summary>
    /// Solves the system of congruences using the Chinese Remainder Theorem (CRT).
    /// </summary>
    /// <param name="a">Array of remainders.</param>
    /// <param name="m">Array of moduli.</param>
    /// <returns>A tuple containing the solution x and the product of all moduli mod.</returns>
    /// <exception cref="ArgumentException">Thrown when the arrays have different lengths or the moduli are not pairwise coprime.</exception>
    public static (BigInteger x, BigInteger mod) CRTbigint(BigInteger[] a, BigInteger[] m)
    {
        if (m.Length < 0 || m.Length != a.Length)
            throw new("array dimension");

        if (m.Any(mi => mi < 2) || m.Grid2D().Where(e => e.t1 != e.t2).Any(e => GcdBigInt(e.t1, e.t2) != 1))
            throw new("modulus");

        var A = a.Zip(m).Select(e => AmodPbigint(e.First, e.Second)).ToArray();
        var mod = m.Aggregate((mi, mj) => mi * mj);
        BigInteger x = 0;
        foreach (var (mi, ai) in m.Zip(A))
        {
            var quo = mod / mi;
            var inv = InvModPbezbigint(quo, mi);
            x = (x + ai * quo * inv) % mod;
        }

        return (x % mod, mod);
    }

    /// <summary>
    /// Solves the equation k^m = 1 (mod n) for a given n and m.
    /// </summary>
    /// <param name="n">The modulus of the equation.</param>
    /// <param name="m">The power of the equation.</param>
    /// <returns>The solution of the equation, or -1 if no solution exists.</returns>
    public static int Solve_k_pow_m_equal_one_mod_n(int n, int m)
    {
        var seq = Enumerable.Range(2, n - 2);
        var criteria = seq.Where(i => Gcd(i, n) == 1 && PowMod(i, m, n) == 1);
        return criteria.FirstOrDefault(-1);
    }

    /// <summary>
    /// Solves the equation k^m = 1 (mod n) for a given n and m.
    /// </summary>
    /// <param name="n">The modulus of the equation.</param>
    /// <param name="m">The power of the equation.</param>
    /// <returns>The solution of the equation, or -1 if no solution exists.</returns>
    public static int Solve_k_pow_m_equal_one_mod_n_strict(int n, int m)
    {
        var seq = Enumerable.Range(1, n - 1).Reverse();
        var criteria = seq.Where(i => Gcd(i, n) == 1 && PowModEqualOne(i, m, n));
        return criteria.FirstOrDefault(-1);
    }

    /// <summary>
    /// Calculates if the result of a^e mod m is equal to 1.
    /// </summary>
    /// <param name="a">The base of the exponent.</param>
    /// <param name="e">The exponent.</param>
    /// <param name="m">The modulus.</param>
    /// <returns>True if a^e mod m is equal to 1, false otherwise.</returns>
    public static bool PowModEqualOne(int a, int e, int m)
    {
        if (PowMod(a, e, m) != 1)
            return false;
        
        var a0 = 1;
        for (var k = 0; k < e; ++k)
        {
            a0 = (a0 * a) % m;
            if (a0 == 1 && k < e - 1)
                return false;
        }

        return a0 == 1;
    }

    /// <summary>
    /// Calculates if the result of a^e mod m is equal to 1.
    /// </summary>
    /// <param name="a">The base of the exponent.</param>
    /// <param name="e">The exponent.</param>
    /// <param name="m">The modulus.</param>
    /// <returns>True if a^e mod m is equal to 1, false otherwise.</returns>
    public static bool PowModEqualOnelong(long a, long e, long m)
    {
        if (PowModLong(a, e, m) != 1)
            return false;

        long a0 = 1;
        for (var k = 0; k < e; ++k)
        {
            a0 = AmodPlong(a0 * a, m);
            if (a0 == 1 && k < e - 1)
                return false;
        }

        return a0 == 1;
    }

    /// <summary>
    /// Calculates if the result of a^e mod m is equal to 1.
    /// </summary>
    /// <param name="a">The base of the exponent.</param>
    /// <param name="e">The exponent.</param>
    /// <param name="m">The modulus.</param>
    /// <returns>True if a^e mod m is equal to 1, false otherwise.</returns>
    public static bool PowModEqualOnebigint(BigInteger a, BigInteger e, BigInteger m)
    {
        if (BigInteger.ModPow(a, e, m) != 1)
            return false;

        BigInteger a0 = 1;
        for (var k = 0; k < e; ++k)
        {
            a0 = AmodPbigint(a0 * a, m);
            if (a0 == 1 && k < e - 1)
                return false;
        }

        return a0 == 1;
    }

    /// <summary>
    /// Calculates the remainder of a divided by p.
    /// </summary>
    /// <param name="a">The dividend.</param>
    /// <param name="p">The divisor.</param>
    /// <returns>The remainder of a divided by p.</returns>
    public static int AmodP(int a, int p)
    {
        int r = a % p;
        return r < 0 ? r + p : r;
    }

    /// <summary>
    /// Calculates the remainder of a divided by p.
    /// </summary>
    /// <param name="a">The dividend.</param>
    /// <param name="p">The divisor.</param>
    /// <returns>The remainder of a divided by p.</returns>
    public static long AmodPlong(long a, long p)
    {
        long r = a % p;
        return r < 0 ? r + p : r;
    }

    /// <summary>
    /// Calculates the remainder of a divided by p.
    /// </summary>
    /// <param name="a">The dividend.</param>
    /// <param name="p">The divisor.</param>
    /// <returns>The remainder of a divided by p.</returns>
    public static BigInteger AmodPbigint(BigInteger a, BigInteger p)
    {
        BigInteger r = a % p;
        return r < 0 ? r + p : r;
    }

    /// <summary>
    /// Calculates the inverse modulo of a number a modulo p.
    /// Using Bezout method
    /// </summary>
    /// <param name="a">The number to calculate the inverse modulo of.</param>
    /// <param name="p">The modulo.</param>
    /// <returns>The inverse modulo of a modulo p.</returns>
    public static int InvModPbez(int a, int p) => AmodP(Bezout(a, p).x, p);

    /// <summary>
    /// Calculates the inverse modulo of a number a modulo p.
    /// Using Bezout method
    /// </summary>
    /// <param name="a">The number to calculate the inverse modulo of.</param>
    /// <param name="p">The modulo.</param>
    /// <returns>The inverse modulo of a modulo p.</returns>
    public static long InvModPbezlong(long a, long p) => AmodPlong(BezoutLong(a, p).x, p);

    /// <summary>
    /// Calculates the inverse modulo of a number a modulo p.
    /// Using Bezout method
    /// </summary>
    /// <param name="a">The number to calculate the inverse modulo of.</param>
    /// <param name="p">The modulo.</param>
    /// <returns>The inverse modulo of a modulo p.</returns>
    public static BigInteger InvModPbezbigint(BigInteger a, BigInteger p) => AmodPbigint(BezoutBigInt(a, p).Xa, p);

    /// <summary>
    /// Creates a dictionary with keys from 1 to n, and values are the invert mod n of the keys. 
    /// </summary>
    /// <param name="n">The upper limit of the range of numbers.</param>
    /// <returns>A Dictionary with keys from 1 to n, and values are the invert mod n of the keys.</returns>
    public static Dictionary<int, int> UnInvertible(int n)
    {
        return (n - 1).SeqLazy(1).Where(i => Gcd(i, n) == 1).ToDictionary(i => i, i => InvModPbez(i, n));
    }

    /// <summary>
    /// Solves the equation k^m = 1 (mod n) and gcd(k, n) = 1.
    /// </summary>
    /// <param name="n">The modulus.</param>
    /// <param name="m">The exponent.</param>
    /// <returns>A list of all solutions to the equation k^m = 1 (mod n).</returns>
    public static IEnumerable<int> SolveAll_k_pow_m_equal_one_mod_n(int n, int m)
    {
        var seq = Enumerable.Range(2, n - 2);
        var criteria = seq.Where(i => Gcd(i, n) == 1 && PowMod(i, m, n) == 1);
        return criteria;
    }

    /// <summary>
    /// Solves the equation k^m = 1 (mod n) and gcd(k, n) = 1.
    /// </summary>
    /// <param name="n">The modulus.</param>
    /// <param name="m">The exponent.</param>
    /// <returns>A list of all solutions to the equation k^m = 1 (mod n).</returns>
    public static IEnumerable<int> SolveAll_k_pow_m_equal_one_mod_n_strict(int n, int m)
    {
        var seq = Enumerable.Range(2, n - 2);
        var criteria = seq.Where(i => Gcd(i, n) == 1 && PowModEqualOne(i, m, n));
        return criteria;
    }

    /// <summary>
    /// Solves the equation k^m = 1 (mod n) and gcd(k, n) = 1.
    /// </summary>
    /// <param name="n">The modulus.</param>
    /// <param name="m">The exponent.</param>
    /// <returns>A list of all solutions to the equation k^m = 1 (mod n).</returns>
    public static IEnumerable<int> SolveAll_k_pow_m_equal_one_mod_n_strict_long(int n, int m)
    {
        var seq = Enumerable.Range(2, n - 2);
        var criteria = seq.Where(i => GcdLong(i, n) == 1 && PowModEqualOnelong(i, m, n));
        return criteria;
    }

    /// <summary>
    /// Solves the equation a^m = 1 (mod n).
    /// </summary>
    /// <param name="n">The modulus.</param>
    /// <param name="a">The base of the equation.</param>
    /// <returns>A sequence of integers representing the solutions to the equation.</returns>
    public static IEnumerable<int> SolveAll_a_pow_m_equal_one_mod_n(int n, int a)
    {
        var seq = Enumerable.Range(2, n - 2);
        return seq.Where(m => PowMod(a, m, n) == 1);
    }

    /// <summary>
    /// Returns the coprimes of a given integer. 
    /// </summary>
    /// <param name="n">The integer to find the coprimes of.</param>
    /// <returns>An IEnumerable containing all the coprimes of n.</returns>
    public static IEnumerable<int> Coprimes(int n) => n.Range(1).Where(a => Gcd(a, n) == 1);

    /// <summary>
    /// Computes the Euler's totient function phi(n) for a given integer n.
    /// </summary>
    /// <param name="n">The integer for which to compute phi(n).</param>
    /// <returns>The Euler's totient function phi(n).</returns>
    public static int Phi(int n) => Coprimes(n).Count();

    /// <summary>
    /// Calculates the dividors of a given integer.
    /// </summary>
    /// <param name="n">The integer to calculate dividors for.</param>
    /// <returns>An enumerable containing the dividors of the given integer.</returns>
    public static IEnumerable<int> Dividors(int n) => Enumerable.Range(1, n / 2).Where(i => i != n && n % i == 0);

    /// <summary>
    /// Calculates the dividors of a given integer.
    /// </summary>
    /// <param name="n">The integer to calculate dividors for.</param>
    /// <returns>An enumerable containing the dividors of the given integer.</returns>
    public static IEnumerable<BigInteger> DividorsBigInt(BigInteger n)
    {
        var decomp = PrimesDecompositionBigInt(BigInteger.Abs(n)).GroupBy(e => e)
            .ToDictionary(e => e.Key, e => e.Count());

        var set = decomp.Select(e => (e.Value + 1).Range().Select(k => BigInteger.Pow(e.Key, k))).MultiLoop();
        return set.Select(l => l.Aggregate(BigInteger.One, (acc, e) => acc * e));
    }

    /// <summary>
    /// Calculates the Carmichael numbers up to the given value. 
    /// </summary>
    /// <param name="n">The maximum value of the Carmichael numbers.</param>
    /// <returns>A list of all Carmichael numbers up to the given value.</returns>
    public static List<int> Carmichael(int n)
    {
        var l = Coprimes(n)
            .Select(a => SolveAll_a_pow_m_equal_one_mod_n(n, a))
            .Aggregate((a, b) => a.Intersect(b));
        return l.ToList();
    }

    /// <summary>
    /// Calculates the least common multiple of two integers.
    /// </summary>
    /// <param name="a">The first integer.</param>
    /// <param name="b">The second integer.</param>
    /// <returns>The least common multiple of the two integers.</returns>
    public static int Lcm(int a, int b)
    {
        return a * b / Gcd(a, b);
    }

    /// <summary>
    /// Calculates the least common multiple of an array of integers.
    /// </summary>
    /// <param name="arr">An array of integers.</param>
    /// <returns>The least common multiple of the integers in the array.</returns>
    public static int Lcm(int[] arr)
    {
        if (!arr.Any())
            return 1;

        return Lcm(arr.First(), Lcm(arr.Skip(1).ToArray()));
    }

    /// <summary>
    /// Generates an array of integers with a range from the given start and step values.
    /// </summary>
    /// <param name="a">The end value of the range.</param>
    /// <param name="start">The starting value of the range (defaults to 0).</param>
    /// <param name="step">The increment between each value (defaults to 1).</param>
    /// <returns>An array of integers with a range from the given start and step values.</returns>
    public static int[] Range(this int a, int start = 0, int step = 1) =>
        Enumerable.Range(0, a).Select(i => start + i * step).ToArray();

    /// <summary>
    /// Calculates and returns the factorial of the given integer.
    /// </summary>
    /// <param name="i">The integer to calculate the factorial of.</param>
    /// <returns>The factorial of the given integer.</returns>
    public static BigInteger Fact(this int i) => Enumerable.Range(0, i).Aggregate(BigInteger.One, (acc, k) => acc * k);

    public static IEnumerable<int> SeqLazy(this int a, int start = 0, int step = 1)
    {
        return Enumerable.Range(0, a).Select(i => start + i * step);
    }

    /// <summary>
    /// Generates all possible permutations of the numbers from 1 to n.
    /// </summary>
    /// <param name="n">The upper bound of the range of numbers.</param>
    /// <returns>An array containing all possible permutations of the numbers from 1 to n.</returns>
    public static int[][] GetPermutations(int n)
    {
        return AllPermutations[n];
    }

    /// <summary>
    /// Gets the kth permutation of a sequence of numbers from 1 to n.
    /// </summary>
    /// <param name="n">The length of the sequence.</param>
    /// <param name="k">The index of the permutation to be returned.</param>
    /// <returns>An array containing the kth permutation of a sequence from 1 to n.</returns>
    public static int[] GetPermutation(int n, int k)
    {
        return AllPermutations[n][k];
    }

    /// <summary>
    /// Generates a hash value for the given integer and array.
    /// </summary>
    /// <param name="n">An integer value.</param>
    /// <param name="m">An array of integers.</param>
    /// <returns>An integer representing the generated hash value.</returns>
    public static int GenHash(int n, int[] m)
    {
        var pow = 1;
        var hash = 0;
        for (var k = 0; k < m.Length; ++k)
        {
            hash += pow * m[k];
            pow *= n;
        }

        return hash;
    }

    /// <summary>
    /// Converts a permutation of length n into its corresponding cycle representation.
    /// </summary>
    /// <param name="n">The length of the permutation.</param>
    /// <param name="p">The permutation to convert.</param>
    /// <returns>An array of cycles representing the given permutation.</returns>
    public static int[][] PermutationToCycles(int n, int[] p)
    {
        var orbits = new List<List<int>>();
        if (p.All(i => i >= 0 && i < n) && n == p.Length && n == p.Distinct().Count())
        {
            var dic = Enumerable.Range(0, n).ToDictionary(i => i, i => p[i]);
            while (dic.Count != 0)
            {
                var v = dic.MinBy(e => e.Key);
                if (v.Key == v.Value)
                {
                    orbits.Add(new List<int> { v.Key });
                    dic.Remove(v.Key);
                    continue;
                }

                var idx = v.Key;
                var cycle = new List<int> { idx };
                while (true)
                {
                    var idx0 = dic[idx];
                    if (idx0 == v.Key)
                    {
                        dic.Remove(idx);
                        break;
                    }

                    cycle.Add(idx0);
                    dic.Remove(idx);
                    idx = idx0;
                }

                orbits.Add(cycle);
            }
        }

        return orbits.Select(l => l.ToArray()).ToArray();
    }

    /// <summary>
    /// Inverts the given permutation of an array.
    /// </summary>
    /// <param name="arr0">The original array.</param>
    /// <param name="arr1">The permutation of the original array.</param>
    /// <returns>An integer representing the inversion of the given permutation.</returns>
    public static int InvertPermutation(int[] arr0, int[] arr1)
    {
        var n = arr0.Length;
        for (var k = 0; k < n; ++k)
            arr1[arr0[k]] = k;

        return GenHash(n, arr1);
    }

    /// <summary>
    /// Compose a permutation from two arrays.
    /// </summary>
    /// <param name="arr0">The first array.</param>
    /// <param name="arr1">The second array.</param>
    /// <param name="arr2">The result array.</param>
    /// <returns>The composed permutation.</returns>
    public static int ComposePermutation(int[] arr0, int[] arr1, int[] arr2)
    {
        var n = arr0.Length;
        var hash = 0;
        var pow = 1;
        for (var k = 0; k < n; ++k)
        {
            var v = arr1[arr0[k]];
            hash += v * pow;
            pow *= n;
            arr2[k] = v;
        }

        return hash;
    }

    /// <summary>
    /// Checks if the given array is a valid permutation of size n.
    /// </summary>
    /// <param name="n">The size of the permutation.</param>
    /// <param name="arr">The array to be checked.</param>
    /// <returns>True if the array is a valid permutation, false otherwise.</returns>
    public static bool CheckTable(int n, int[] arr)
    {
        if (arr.Length != n)
            return false;

        if (arr.Min() < 0 || arr.Max() > n - 1 || arr.Distinct().Count() != n)
            return false;

        return true;
    }

    /// <summary>
    /// Checks if a cycle exists in the given array.
    /// </summary>
    /// <param name="n">Number of elements in the array.</param>
    /// <param name="arr">Array to be checked for a cycle.</param>
    /// <returns>True if a cycle exists, false otherwise.</returns>
    public static bool CheckCycle(int n, int[] arr)
    {
        if (arr.Min() < 0 || arr.Max() > n - 1 || arr.Distinct().Count() != arr.Length)
            return false;

        return true;
    }

    /// <summary>
    /// Applies the given cycle to the given array.
    /// </summary>
    /// <param name="arr">The array to apply the cycle to.</param>
    /// <param name="cycle">The cycle to apply.</param>
    public static void ApplyCycle(int[] arr, int[] cycle)
    {
        var a0 = arr[cycle[0]];
        var n = cycle.Length;
        for (var k = 0; k < n - 1; ++k)
            arr[cycle[k]] = arr[cycle[k + 1]];

        arr[cycle[n - 1]] = a0;
    }

    /// <summary>
    /// Gets a permutation and its cycles representation from a given type
    /// </summary>
    /// <param name="type">The type of the permutation.</param>
    /// <returns>A permutation and an array of cycles.</returns>
    public static (int[] perm, int[][] cycles) PermAndCyclesFromType(int[] type)
    {
        var dim = type.Sum();
        var lt = new List<int[]>();
        var rg = dim.Range();
        var perm = new int[dim];
        foreach (var k in type.Order())
        {
            var r0 = rg.Take(k).ToArray();
            rg = rg.Skip(k).ToArray();
            lt.Add(r0);
            for (int i = 0; i < k; i++)
                perm[r0[i]] = r0[(i + 1) % k];
        }

        return (perm, lt.ToArray());
    }

    /// <summary>
    /// Calculates the nth root of a given number.
    /// </summary>
    /// <param name="y">The number for which the nth root is to be calculated.</param>
    /// <param name="n">The degree of the root to be calculated.</param>
    /// <returns>The nth root of the given number.</returns>
    public static BigInteger NthRootBigint(BigInteger y, int n)
    {
        if (n <= 0)
            return 1;

        if (n % 2 == 0 && y.Sign == -1)
            throw new("Even NthRoot must has positive argument");

        if (n == 1)
            return y;

        if (y.IsOne)
            return y;

        BigInteger ai;
        var approx = 10;
        var ya = y * BigInteger.Pow(approx, n);
        var aj = 1 + ya / n;
        do
        {
            ai = aj;
            var aiPow = BigInteger.Pow(ai, n - 1);
            var num = aiPow * ai - ya;
            var denom = n * aiPow;
            aj = ai - num / denom; // Newton iteration
        } while (BigInteger.Abs(ai - aj) > approx / 3);

        return aj / approx;
    }

    /// <summary>
    /// Calculates the square root of a given number.
    /// </summary>
    /// <param name="y">The number for which the square root is to be calculated.</param>
    /// <returns>The square root of the given number.</returns>
    public static BigInteger SqrtBigInt(BigInteger y) => NthRootBigint(y, 2);
}