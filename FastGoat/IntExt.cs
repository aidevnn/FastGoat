namespace FastGoat;

public static class IntExt
{
    static IntExt()
    {
        AllPermutations = new();
        var acc = new List<List<int>>() { new List<int>() };
        var comp = Comparer<List<int>>.Create((a, b) => a.SequenceCompare(b));
        for (int i = 0; i < 8; ++i)
        {
            var tmpAcc = new List<List<int>>();
            foreach (var l0 in acc)
            {
                for (int k = 0; k <= i; ++k)
                {
                    var l1 = l0.ToList();
                    l1.Insert(k, i);
                    tmpAcc.Add(l1);
                }
            }

            acc = tmpAcc.ToList();
            AllPermutations[i + 1] = acc.OrderBy(a => a, comp).Select(a => a.Select(b => b + 1).ToArray()).ToArray();
        }
    }
    static Dictionary<int, int[][]> AllPermutations { get; }
    public static int GCD(int a, int b)
    {
        if (a < b)
            return GCD(b, a);

        if (b == 0)
            return a;

        int q = a / b;
        int r = a - b * q;
        return GCD(b, r);
    }
    public static int PowMod(int a, int exp, int mod)
    {
        var a0 = 1;
        for (int k = 0; k < exp; ++k)
            a0 = (a0 * a) % mod;

        return a0;
    }

    public static int LCM(int a, int b) => a * b / GCD(a, b);
    public static int[] Range(this int a, int start = 0) => Enumerable.Range(start, a).ToArray();
    public static int[][] GetPermutations(int n) => AllPermutations[n];
    public static int[] GetPermutation(int n, int k) => AllPermutations[n][k];
    public static int GenHash(int n, int[] m)
    {
        var pow = 1;
        var hash = 0;
        for (int k = 0; k < m.Length; ++k)
        {
            hash += pow * m[k];
            pow *= n;
        }

        return hash;
    }
    public static HashSet<HashSet<int>> Orbits(this int[] arr)
    {
        HashSet<HashSet<int>> hs = new HashSet<HashSet<int>>(new EqualityHashSet<int>());
        for (int k = 0; k < arr.Length; ++k)
        {
            var a0 = k;
            int sz = 0;
            var hs0 = new HashSet<int>() { a0 };
            while (sz != hs0.Count)
            {
                sz = hs0.Count;
                a0 = arr[a0];
                hs0.Add(a0);
            }
            hs.Add(hs0);
        }

        return hs;
    }

    public static int InvertPermutation(int[] arr0, int[] arr1)
    {
        int n = arr0.Length;
        for (int k = 0; k < n; ++k)
            arr1[arr0[k]] = k;

        return GenHash(n, arr1);
    }

    public static int ComposePermutation(int[] arr0, int[] arr1, int[]? arr2 = null)
    {
        int n = arr0.Length;
        int hash = 0;
        int pow = 1;
        for (int k = 0; k < n; ++k)
        {
            var v = arr1[arr0[k]];
            hash += v * pow;
            pow *= n;
            arr2?.SetValue(v, k);
        }

        return hash;
    }

    public static void Add(this int[] arr, int v)
    {
        for (int k = 0; k < arr.Length; ++k)
            arr[k] += v;
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
        int a0 = arr[cycle[0]];
        int n = cycle.Length;
        for (int k = 0; k < n - 1; ++k)
            arr[cycle[k]] = arr[cycle[k + 1]];

        arr[cycle[n - 1]] = a0;
    }

}
