namespace FastGoat;

public static class IntExt
{
    static IntExt()
    {
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
            AllPermutations[i + 1] = acc.OrderBy(a => a, comp).Select(a => a.Select(b => b + 1).ToArray()).ToArray();
        }
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

    public static int Solve_x_pow_n_equal_one_mod_m(int m, int n)
    {
        var seq = Enumerable.Range(2, m - 2);
        var criteria = seq.Where(i => Gcd(i, m) == 1 && PowMod(i, n, m) == 1);
        return criteria.FirstOrDefault();
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
                    if (idx0 == v.Value)
                    {
                        dic.Remove(v.Key);
                        break;
                    }

                    cycle.Add(idx0);
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