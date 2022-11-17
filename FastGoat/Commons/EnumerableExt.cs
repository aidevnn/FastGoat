namespace FastGoat.Commons;

public static class EnumerableExt
{
    public static string Glue<T>(this IEnumerable<T> ts, string sep = "", string fmt = "{0}")
    {
        return string.Join(sep, ts.Select(t => string.Format(fmt, t)));
    }
    //
    // public static string GlueMap<T1, T2>(this IDictionary<T1, T2> map, string sep = ", ", string fmt = "{0}->{1}")
    // {
    //     return string.Join(sep, map.Select(kp => string.Format(fmt, kp.Key, kp.Value)));
    // }

    public static string GlueMap<T1, T2>(this IEnumerable<KeyValuePair<T1, T2>> map, string sep = ", ",
        string fmt = "{0}->{1}")
    {
        return string.Join(sep, map.Select(kp => string.Format(fmt, kp.Key, kp.Value)));
    }

    public static IOrderedEnumerable<T> Ascending<T>(this IEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.OrderBy(t => t);
    }

    public static IOrderedEnumerable<IEnumerable<T>> AscendingByCount<T>(this IEnumerable<IEnumerable<T>> ts)
    {
        return ts.OrderBy(t => t.Count());
    }

    public static IOrderedEnumerable<KeyValuePair<T1, T2>> AscendingByKey<T1, T2>(
        this IEnumerable<KeyValuePair<T1, T2>> ts) where T1 : IComparable<T1>
    {
        return ts.OrderBy(t => t.Key);
    }

    public static IOrderedEnumerable<T> ThenAscending<T>(this IOrderedEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.ThenBy(t => t);
    }

    public static IOrderedEnumerable<KeyValuePair<T1, T2>> ThenAscendingByKey<T1, T2>(
        this IOrderedEnumerable<KeyValuePair<T1, T2>> ts) where T1 : IComparable<T1>
    {
        return ts.ThenBy(t => t.Key);
    }

    public static IOrderedEnumerable<IEnumerable<T>> ThenAscendingByCount<T>(this IOrderedEnumerable<IEnumerable<T>> ts)
    {
        return ts.ThenBy(t => t.Count());
    }

    public static IOrderedEnumerable<T> Descending<T>(this IEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.OrderByDescending(t => t);
    }

    public static IOrderedEnumerable<KeyValuePair<T1, T2>> DescendingByKey<T1, T2>(
        this IEnumerable<KeyValuePair<T1, T2>> ts) where T1 : IComparable<T1>
    {
        return ts.OrderByDescending(t => t.Key);
    }

    public static IOrderedEnumerable<IEnumerable<T>> DescendingByCount<T>(this IEnumerable<IEnumerable<T>> ts)
    {
        return ts.OrderByDescending(t => t.Count());
    }

    public static IOrderedEnumerable<T> ThenDescending<T>(this IOrderedEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.ThenByDescending(t => t);
    }

    public static IOrderedEnumerable<KeyValuePair<T1, T2>> ThenDescendingByKey<T1, T2>(
        this IOrderedEnumerable<KeyValuePair<T1, T2>> ts) where T1 : IComparable<T1>
    {
        return ts.ThenByDescending(t => t.Key);
    }

    public static IOrderedEnumerable<IEnumerable<T>> ThenDescendingByCount<T>(
        this IOrderedEnumerable<IEnumerable<T>> ts)
    {
        return ts.ThenByDescending(t => t.Count());
    }

    public static int SequenceCompareTo<T>(this IEnumerable<T> sq1, IEnumerable<T> sq2) where T : IComparable<T>
    {
        var en1 = sq1.GetEnumerator();
        var en2 = sq2.GetEnumerator();
        while (en1.MoveNext() && en2.MoveNext())
        {
            var e1 = en1.Current;
            var e2 = en2.Current;
            var comp1 = e1.CompareTo(e2);
            if (comp1 == 0)
                continue;

            en1.Dispose();
            en2.Dispose();
            return comp1;
        }

        en1.Dispose();
        en2.Dispose();
        // return sq1.Count().CompareTo(sq2.Count());
        return 0;
    }

    // x3 Faster than append concat, for permutation methods
    public static IEnumerable<T> InsertAt<T>(this IEnumerable<T> ts, int index, T p)
    {
        int k = 0;
        foreach (var t in ts)
        {
            if (index == k++)
                yield return p;

            yield return t;
        }

        if (index == k)
            yield return p;
    }

    public static IEnumerable<IEnumerable<T>> AllCombinations<T>(this T[] seq)
    {
        if (seq.Length > IntExt.NbCombinations)
            throw new Exception($"Max length is {IntExt.NbCombinations}");

        return IntExt.GetCombinations(seq.Length)
            .Select(comb => comb.Zip(seq).Where(e => e.First).Select(e => e.Second));
    }

    public static IEnumerable<(T1, T2)> MultiLoop<T1, T2>(IEnumerable<T1> t1s, IEnumerable<T2> t2s)
    {
        return t1s.SelectMany(t1 => t2s.Select(t2 => (t1, t2)));
    }

    public static IEnumerable<IEnumerable<T>> MultiLoop<T>(this IEnumerable<IEnumerable<T>> enumerables)
    {
        if (!enumerables.Any())
            yield return Enumerable.Empty<T>(); // on the void
        else
        {
            foreach (var t in enumerables.First())
            {
                foreach (var v in MultiLoop<T>(enumerables.Skip(1)))
                {
                    yield return v.Prepend(t);
                }
            }
        }
    }

    public static IEnumerable<IEnumerable<T>> MultiLoopWith<T>(this IEnumerable<T> ts, params IEnumerable<T>[] other)
    {
        return MultiLoop(other.Prepend(ts));
    }

    public static IEnumerable<IEnumerable<T>> MultiLoop<T>(this IEnumerable<T> arr, int n)
    {
        return Enumerable.Repeat(arr, n).MultiLoop();
    }

    public static IEnumerable<(T1 t1, T2 t2)> Grid2D<T1, T2>(this IEnumerable<T1> first, IEnumerable<T2> second)
    {
        foreach (var t1 in first)
        {
            foreach (var t2 in second)
            {
                yield return (t1, t2);
            }
        }
    }

    public static IEnumerable<(T t1, T t2)> Grid2D<T>(this T[] seq) => seq.Grid2D(seq);
}

public class SetEquality<T> : EqualityComparer<HashSet<T>> where T : IEquatable<T>
{
    public override bool Equals(HashSet<T>? x, HashSet<T>? y)
    {
        return y is not null && x is not null && x.SetEquals(y);
    }

    public override int GetHashCode(HashSet<T> obj) => (obj.Count(), typeof(T).GetHashCode()).GetHashCode();
}

public class SequenceEquality<T> : EqualityComparer<IEnumerable<T>> where T : IEquatable<T>, IComparable<T>
{
    public override bool Equals(IEnumerable<T>? x, IEnumerable<T>? y)
    {
        return y is not null && x is not null && x.SequenceEqual(y);
    }

    public override int GetHashCode(IEnumerable<T> obj) => (obj.Count(), typeof(T).GetHashCode()).GetHashCode();
}