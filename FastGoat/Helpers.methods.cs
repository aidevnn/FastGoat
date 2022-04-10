using System;
using System.Collections.Generic;
using System.Linq;

namespace FastGoat
{
    public static class Helpers
    {
        public static string Glue<T>(this IEnumerable<T> ts, string fmt = "{0,2}", string sep = " ")
        {
            return string.Join(sep, ts.Select(d => string.Format(fmt, d)));
        }

        public static int[] Range(this int a, int start = 0) => Enumerable.Range(start, a).ToArray();

        public static int Prod(this IEnumerable<int> vs) => vs.Aggregate((a, acc) => a * acc);

        public static int Factoriel(this int n) => n.Range(1).Prod();

        public static int FactorielN(this int n) => n > 12 ? 1000000 : n.Factoriel();

        public static void ClearArray(int[] arr)
        {
            for (int k = 0; k < arr.Length; ++k)
                arr[k] = 0;
        }

        public static int AddModulo(int[] dims, int[] m0, int[] m1, int[] m2 = null)
        {
            int pow = 1;
            int hash = 0;
            for (int k = 0; k < Math.Min(dims.Length, m0.Length); ++k)
            {
                var v = (m0[k] + m1[k]) % dims[k];
                hash += pow * v;
                pow *= dims[k];
                m2?.SetValue(v, k);
            }

            return hash;
        }

        public static int InvertModulo(int[] dims, int[] m0, int[] m1 = null)
        {
            int pow = 1;
            int hash = 0;
            for (int k = 0; k < Math.Min(dims.Length, m0.Length); ++k)
            {
                var n = dims[k];
                var v = (n - m0[k]) % n;
                hash += pow * v;
                pow *= n;
                m1?.SetValue(v, k);
            }

            return hash;
        }

        public static int ArrayCompare(int[] x, int[] y)
        {
            for (int k = 0; k < Math.Min(x.Length, y.Length); ++k)
            {
                var e0 = x[k];
                var e1 = y[k];
                if (e0 != e1)
                    return e0.CompareTo(e1);
            }
            return 0;
        }

        public static int GenHash(int[] n, int[] m)
        {
            var pow = 1;
            var hash = 0;
            for (int k = 0; k < m.Length; ++k)
            {
                hash += pow * m[k];
                pow *= n[k];
            }

            return hash;
        }

        public static int[] PrimeDecomposition(int n)
        {
            var l = new List<int>();
            var n0 = (int)Math.Sqrt(n) + 1;
            int m = n;
            for(int p = 1; p < n0; ++p)
            {
                while (m % p == 0)
                {
                    l.Add(p);
                    m /= p;
                }
            }

            return l.ToArray();
        }

        public static int[] Canonic(int n, int k)
        {
            var arr = new int[n];
            arr[k] = 1;
            return arr;
        }

        public static int[][] BaseCanonic(int n) => n.Range().Select(k => Canonic(n, k)).ToArray();

        public static int[][] AllTuples(params int[] dims)
        {
            var acc = new List<List<int>>() { new List<int>() };
            for (int i = 0; i < dims.Length; ++i)
            {
                var tmpAcc = new List<List<int>>();
                foreach (var l0 in acc)
                {
                    for (int k = 0; k < dims[i]; ++k)
                    {
                        var l1 = l0.ToList();
                        l1.Add(k);
                        tmpAcc.Add(l1);
                    }
                }

                acc = tmpAcc.ToList();
            }

            return acc.Select(a => a.ToArray()).ToArray();
        }

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

        public static int InvertPermutation(int[] arr0, int[] arr1)
        {
            int n = arr0.Length;
            for (int k = 0; k < n; ++k)
                arr1[arr0[k]] = k;

            return GenHash(n, arr1);
        }

        public static int ComposePermutation(int[] arr0, int[] arr1, int[] arr2 = null)
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

        public static int[][] AllPerms(int n)
        {
            var acc = new List<List<int>>() { new List<int>() };
            for (int i = 0; i < n; ++i)
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
            }

            return acc.Select(a => a.ToArray()).ToArray();
        }

        public static HashSet<HashSet<int>> Orbits(this int[] arr)
        {
            HashSet<HashSet<int>> hs = new HashSet<HashSet<int>>(new EqualityHashSet());
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
    }

    public class EqualityHashSet : EqualityComparer<HashSet<int>>
    {
        public override bool Equals(HashSet<int> x, HashSet<int> y) => x.SetEquals(y);
        public override int GetHashCode(HashSet<int> obj) => obj.Count;
    }

}
