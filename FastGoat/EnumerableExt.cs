namespace FastGoat;

public static class EnumerableExt
{
    public static string Glue<U>(this IEnumerable<U> us, string sep = "", string fmt = "{0}")
    {
        return string.Join(sep, us.Select(u => string.Format(fmt, u)));
    }
    public static string GlueMaxDigits<U>(this IEnumerable<U> us, string sep = "")
    {
        var digits = us.Max(u => u?.ToString()?.Length);
        var fmt = $"{{0,{digits}}}";
        return string.Join(sep, us.Select(u => string.Format(fmt, u)));
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
    public static int SequenceCompare<U>(this IEnumerable<U> first, IEnumerable<U> second) where U : IComparable<U>
    {
        var firstEnum = first.GetEnumerator();
        var secondEnum = second.GetEnumerator();
        while (firstEnum.MoveNext() && secondEnum.MoveNext())
        {
            var comp = firstEnum.Current.CompareTo(secondEnum.Current);
            if (comp != 0)
                return comp;
        }

        return first.Count().CompareTo(second.Count());
    }

    public static IEnumerable<(U, V)> Zip2<U, V>(this IEnumerable<U> first, IEnumerable<V> second)
    {
        var firstEnum = first.GetEnumerator();
        var secondEnum = second.GetEnumerator();
        while (firstEnum.MoveNext() && secondEnum.MoveNext())
        {
            yield return (firstEnum.Current, secondEnum.Current);
        }
    }
    public static IEnumerable<Z> Zip2<U, V, Z>(this IEnumerable<U> first, IEnumerable<V> second, Func<U, V, Z> func)
    {
        var firstEnum = first.GetEnumerator();
        var secondEnum = second.GetEnumerator();
        while (firstEnum.MoveNext() && secondEnum.MoveNext())
        {
            yield return func(firstEnum.Current, secondEnum.Current);
        }
    }
}
