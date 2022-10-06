namespace FastGoat
{
    public static class IntExt
    {
        static IntExt()
        {
            Primes10000 = new(AllPrimes(10000));
            Partitions32 = IntPartitions(32);

            AllPermutations = new Dictionary<int, int[][]>();
            var acc = new List<List<int>> { new() };
            var comp = Comparer<List<int>>.Create((a, b) => a.SequenceCompareTo(b));
            for (var i = 0; i < 8; ++i)
            {
                var tmpAcc = new List<List<int>>();
                foreach (var l0 in acc)
                    for (var k = 0; k <= i; ++k)
                    {
                        var l1 = l0.ToList();
                        l1.Insert(k, i);
                        tmpAcc.Add(l1);
                    }

                acc = tmpAcc.ToList();
                AllPermutations[i + 1] =
                    acc.OrderBy(a => a, comp).Select(a => a.Select(b => b + 1).ToArray()).ToArray();
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
        static List<(int p, int pow)> PrimesDecomposition(int n)
        {
            List<(int p, int pow)> decomp = new();
            var n0 = n;
            foreach (var p in Primes10000.Where(e => e <= n))
            {
                var pow = 0;
                while (n0 % p == 0)
                {
                    ++pow;
                    n0 /= p;
                }

                if (pow != 0)
                    decomp.Add((p, pow));

                if (n0 == 1) break;
            }

            return decomp;
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

        private static Dictionary<int, int[][]> AllPermutations { get; }

        public static int Gcd(int a, int b)
        {
            if (a < 0 && b < 0)
                return Gcd(-a, -b);
            if (a * b < 0)
                if (a < 0)
                    return Gcd(-a, b);
                else
                    return Gcd(a, -b);

            if (a < b)
                return Gcd(b, a);

            if (b == 0)
                return a;

            var q = a / b;
            var r = a - b * q;
            return Gcd(b, r);
        }

        public static int PowMod(int a, int exp, int mod)
        {
            var a0 = 1;
            for (var k = 0; k < exp; ++k)
                a0 = a0 * a % mod;

            return a0;
        }

        // Saunders MacLane, Garrett Birkhoff. Algebra (3rd ed.) criteria
        public static int Solve_k_pow_m_equal_one_mod_n(int n, int m)
        {
            var seq = Enumerable.Range(2, n - 2);
            var criteria = seq.Where(i => Gcd(i, n) == 1 && PowMod(i, m, n) == 1);
            return criteria.FirstOrDefault();
        }

        // Saunders MacLane, Garrett Birkhoff. Algebra (3rd ed.) criteria
        public static List<int> SolveAll_k_pow_m_equal_one_mod_n(int n, int m)
        {
            var seq = Enumerable.Range(2, n - 2);
            var criteria = seq.Where(i => Gcd(i, n) == 1 && PowMod(i, m, n) == 1);
            return criteria.ToList();
        }

        public static int Lcm(int a, int b)
        {
            return a * b / Gcd(a, b);
        }

        public static List<(int a, int b, int q, int r, int x, int y)> BezoutDetails(int a, int b)
        {
            List<(int a, int b, int q, int r, int x, int y)> seq = new();
            if (a < 1 || b < 1)
                return seq;

            int x = 0, y = 1;
            while (b != 0)
            {
                // gcd = ax+by
                // gcd = bx+ry
                // gcd = bx+(a-bq)y
                // gcd = ay + b(x-qy)

                var q = a / b;
                var r = a - b * q;
                seq.Add((a, b, q, r, y, x - q * y));
                a = b;
                b = r;
                var tmp = y;
                y = x - q * y;
                x = tmp;
            }

            return seq;
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