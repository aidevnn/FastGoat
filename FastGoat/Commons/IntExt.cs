namespace FastGoat.Commons
{
    public static class IntExt
    {
        static IntExt()
        {
            Primes10000 = new(AllPrimes(10000));
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

            AllCombinations = new Dictionary<int, bool[][]>
            {
                [0] = new[] { Array.Empty<bool>() }
            };
            for (int k = 1; k <= NbCombinations; ++k)
            {
                AllCombinations[k] = Enumerable.Range(0, 2)
                    .SelectMany(c =>
                        AllCombinations[k - 1].Select(l => l.Prepend(c == 1).ToArray()))
                    .ToArray();
            }
        }

        // Primes Sequence
        public static List<int> Primes10000 { get; }

        // Partitions of integers, 
        // for each key, the dictionary contains all partitions
        public static Dictionary<int, List<List<int>>> Partitions32;

        // Compute all primes less than n
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

        // Primes decompositions of an integer N
        // Return a sequence of tuple (p, pow)
        // N = (p0^pow0) x (p1^pow1) x (p2^pow3) x ... x (p_i^pow_i)
        // Example : The result for N=60 is {(2,2), (3,1), (5,1)}
        public static IEnumerable<int> PrimesDecomposition(int n)
        {
            var n0 = n;
            foreach (var p in Primes10000.Where(e => e <= n))
            {
                while (n0 % p == 0)
                {
                    yield return p;
                    n0 /= p;
                }

                if (n0 == 1) break;
            }
        }

        // Partitions of an integers until N
        // Return a dictionary where each key is mapped to a collection
        // of ordered sequence of integers
        // Example : 5 => {5}, {4,1}, {3,2}, {3,1,1}, {2,2,1}, {2,1,1,1}, {1,1,1,1,1}
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

        public const int NbPermutations = 8;
        public const int NbCombinations = 10;
        private static Dictionary<int, int[][]> AllPermutations { get; }
        private static Dictionary<int, bool[][]> AllCombinations { get; }

        public static IEnumerable<IEnumerable<int>> YieldAllPermutations(int n)
        {
            if (n == 0)
                yield return Enumerable.Empty<int>();
            else
                foreach (var perm in YieldAllPermutations(n - 1))
                    for (int i = 0; i < n; i++)
                        yield return perm.InsertAt(i, n);
        }

        public static IEnumerable<IEnumerable<bool>> YieldAllCombinations(int n)
        {
            return new[] { false, true }.MultiLoop(n);
        }

        public static int Gcd(int a, int b)
        {
            if (b == 0)
                return a < 0 ? -a : a;

            var q = a / b;
            var r = a - b * q;
            return Gcd(b, r);
        }

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

        public static int PowMod(int a, int exp, int mod)
        {
            var a0 = 1;
            for (var k = 0; k < exp; ++k)
                a0 = a0 * a % mod;

            return a0;
        }

        public static IEnumerable<int> LoopPowMod(int a, int mod)
        {
            HashSet<int> set = new() { a };
            int a0 = a;
            while (true)
            {
                a0 = (a0 * a) % mod;
                if (!set.Add(a0))
                    break;
            }

            return set;
        }

        // Saunders MacLane, Garrett Birkhoff. Algebra (3rd ed.) criteria
        public static int Solve_k_pow_m_equal_one_mod_n(int n, int m)
        {
            var seq = Enumerable.Range(2, n - 2);
            var criteria = seq.Where(i => Gcd(i, n) == 1 && PowMod(i, m, n) == 1);
            return criteria.FirstOrDefault();
        }

        public static int KpowMmodN(int n, int m)
        {
            var seq = Enumerable.Range(2, n - 2);
            var criteria = seq.Where(i => Gcd(i, n) == 1 && PowMod(i, m, n) == 1);
            return criteria.FirstOrDefault();
        }

        public static int AmodP(int a, int p)
        {
            int r = a % p;
            return r < 0 ? r + p : r;
        }

        public static int InvModP(int a, int p) => Enumerable.Range(1, p - 1).First(e => AmodP(e * a, p) == 1);

        public static Dictionary<int, int> UnInvertible(int n)
        {
            var seq = from a in Enumerable.Range(1, n)
                from b in Enumerable.Range(1, n)
                select (a, b);
            return seq.Where(i => (i.a * i.b) % n == 1).ToDictionary(i => i.a, i => i.b);
        }

        // Saunders MacLane, Garrett Birkhoff. Algebra (3rd ed.) criteria
        public static List<int> SolveAll_k_pow_m_equal_one_mod_n(int n, int m)
        {
            var seq = Enumerable.Range(2, n - 2);
            var criteria = seq.Where(i => Gcd(i, n) == 1 && PowMod(i, m, n) == 1);
            return criteria.ToList();
        }

        public static IEnumerable<int> SolveAll_a_pow_m_equal_one_mod_n(int n, int a)
        {
            var seq = Enumerable.Range(2, n - 2);
            return seq.Where(m => PowMod(a, m, n) == 1);
        }

        public static IEnumerable<int> Coprimes(int n) => n.Range(1).Where(a => Gcd(a, n) == 1);
        public static int Phi(int n) => Coprimes(n).Count();

        public static List<int> Carmichael(int n)
        {
            var l = Coprimes(n)
                .Select(a => SolveAll_a_pow_m_equal_one_mod_n(n, a))
                .Aggregate((a, b) => a.Intersect(b));
            return l.ToList();
        }

        public static int Lcm(int a, int b)
        {
            return a * b / Gcd(a, b);
        }

        public static int[] Range(this int a, int start = 0)
        {
            return Enumerable.Range(start, a).ToArray();
        }

        public static int[][] GetPermutations(int n)
        {
            return AllPermutations[n];
        }

        public static int[] GetPermutation(int n, int k)
        {
            return AllPermutations[n][k];
        }

        public static bool[][] GetCombinations(int n)
        {
            return AllCombinations[n];
        }

        public static bool[] GetCombination(int n, int k)
        {
            return AllCombinations[n][k];
        }

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

        public static int InvertPermutation(int[] arr0, int[] arr1)
        {
            var n = arr0.Length;
            for (var k = 0; k < n; ++k)
                arr1[arr0[k]] = k;

            return GenHash(n, arr1);
        }

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

        public static bool CheckTable(int n, int[] arr)
        {
            if (arr.Length != n)
                return false;

            if (arr.Min() < 0 || arr.Max() > n - 1 || arr.Distinct().Count() != n)
                return false;

            return true;
        }

        public static bool CheckCycle(int n, int[] arr)
        {
            if (arr.Min() < 0 || arr.Max() > n - 1 || arr.Distinct().Count() != arr.Length)
                return false;

            return true;
        }

        public static void ApplyCycle(int[] arr, int[] cycle)
        {
            var a0 = arr[cycle[0]];
            var n = cycle.Length;
            for (var k = 0; k < n - 1; ++k)
                arr[cycle[k]] = arr[cycle[k + 1]];

            arr[cycle[n - 1]] = a0;
        }
    }
}