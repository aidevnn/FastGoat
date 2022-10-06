namespace FastGoat;

public static class EnumerableExt
{
    public static string Glue<T>(this IEnumerable<T> ts, string sep = "", string fmt = "{0}")
    {
        return string.Join(sep, ts.Select(t => string.Format(fmt, t)));
    }

    public static string GlueMap<T1, T2>(this IDictionary<T1, T2> map, string sep = ", ", string fmt = "{0}->{1}")
    {
        return string.Join(sep, map.Select(kp => string.Format(fmt, kp.Key, kp.Value)));
    }
    public static IOrderedEnumerable<T> Ascending<T>(this IEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.OrderBy(t => t);
    }

    public static IOrderedEnumerable<T> ThenAscending<T>(this IOrderedEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.ThenBy(t => t);
    }

    public static IOrderedEnumerable<T> Descending<T>(this IEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.OrderByDescending(t => t);
    }

    public static IOrderedEnumerable<T> ThenDescending<T>(this IOrderedEnumerable<T> ts) where T : IComparable<T>
    {
        return ts.ThenByDescending(t => t);
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
}