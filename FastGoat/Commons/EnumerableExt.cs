using System.Collections;
using System.Runtime.CompilerServices;

namespace FastGoat.Commons;

public static class EnumerableExt
{
    /// <summary>
    /// Join into a string the elements of an <see cref="IEnumerable{T}"/> together with a separator and format string.
    /// </summary>
    /// <typeparam name="T">The type of elements in the enumerable.</typeparam>
    /// <param name="ts">The enumerable to glue.</param>
    /// <param name="sep">The separator to use between elements. Default is an empty string.</param>
    /// <param name="fmt">The format string to use for each element. Default is "{0}".</param>
    /// <returns>A string containing the joined elements.</returns>
    public static string Glue<T>(this IEnumerable<T> ts, string sep = "", string fmt = "{0}")
    {
        return string.Join(sep, ts.Select(t => string.Format(fmt, t)));
    }

    /// <summary>
    /// Prints the object to the console with a header and a format. 
    /// </summary>
    /// <param name="v">The object to be printed.</param>
    /// <param name="header">The header that will be printed before the object.</param>
    /// <param name="fmt">The format of the output.</param>
    public static void Println(this object v, string header = "Lines", string fmt = "    {0}")
    {
        Console.WriteLine(header);
        if (v is ITuple v0)
            Console.WriteLine(v0.Length.Range().Select(i => v0[i]).Glue("\n", fmt));
        else if (v is IEnumerable v1)
            Console.WriteLine(v1.Cast<object>().Glue("\n", fmt));
        else
            Console.WriteLine(fmt, v);
    }

    /// <summary>
    /// Join into a string a map of type <typeparamref name="T1"/> and <typeparamref name="T2"/>.
    /// </summary>
    /// <typeparam name="T1">The type of the keys in the map.</typeparam>
    /// <typeparam name="T2">The type of the values in the map.</typeparam>
    /// <param name="map">The map to be glued.</param>
    /// <param name="sep">The separator between each key-value pair, defaults to ", ".</param>
    /// <param name="fmt">The format of each key-value pair, defaults to "{0}->{1}".</param>
    /// <returns>A string representation of the given map.</returns>
    public static string GlueMap<T1, T2>(this IEnumerable<KeyValuePair<T1, T2>> map, string sep = ", ", string fmt = "{0}->{1}")
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

    /// <summary>
    /// Compares two sequences of elements of type T which implement the IComparable<T> interface.
    /// </summary>
    /// <typeparam name="T">The type of elements in the sequences.</typeparam>
    /// <param name="sq1">The first sequence to compare.</param>
    /// <param name="sq2">The second sequence to compare.</param>
    /// <returns>A positive value if sq1 is greater than sq2, a negative value if sq1 is less than sq2, or 0 if they are equal.</returns>
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

    /// <summary>
    /// Inserts an element at a specified index in an IEnumerable.
    /// </summary>
    /// <typeparam name="T">The type of the elements in the IEnumerable.</typeparam>
    /// <param name="ts">The IEnumerable to insert the element into.</param>
    /// <param name="index">The index at which to insert the element.</param>
    /// <param name="p">The element to insert.</param>
    /// <returns>An IEnumerable containing the inserted element.</returns>
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

    /// <summary>
    /// Generates all possible combinations of elements in a given sequence.
    /// </summary>
    /// <typeparam name="T">The type of the elements in the sequence.</typeparam>
    /// <param name="seq">The sequence of elements to generate combinations from.</param>
    /// <returns>An enumerable of arrays containing all possible combinations of elements from the given sequence.</returns>
    public static IEnumerable<T[]> AllCombinations<T>(this IEnumerable<T> seq)
    {
        var enumerable = seq as T[] ?? seq.ToArray();
        var nb = enumerable.Count();
        return IntExt.YieldAllCombs(nb).Select(comb => comb.Zip(enumerable).Where(e => e.First).Select(e => e.Second).ToArray());
    }

    /// <summary>
    /// Generates a sequence of sequences by combining elements from the given enumerables.
    /// </summary>
    /// <typeparam name="T">The type of the elements in the enumerables.</typeparam>
    /// <param name="enumerables">The enumerables to combine.</param>
    /// <returns>A sequence of sequences, each containing one element from each of the given enumerables.</returns>
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

    /// <summary>
    /// Generates a sequence of sequences by combining elements from the given IEnumerable and an array of other IEnumerables.
    /// </summary>
    /// <typeparam name="T">The type of elements in the IEnumerable.</typeparam>
    /// <param name="ts">The original IEnumerable.</param>
    /// <param name="other">An array of other IEnumerables.</param>
    /// <returns>A sequence of sequences.</returns>
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

    public static IEnumerable<(T1 t1, T2 t2, T3 t3)> Grid3D<T1, T2, T3>(this IEnumerable<T1> first, IEnumerable<T2> second,
        IEnumerable<T3> third)
    {
        foreach (var t1 in first)
        {
            foreach (var t2 in second)
            {
                foreach (var t3 in third)
                {
                    yield return (t1, t2, t3);
                }
            }
        }
    }

    public static IEnumerable<(T t1, T t2, T t3)> Grid3D<T>(this T[] seq) => Grid3D(seq, seq, seq);
}

/// <summary>
/// Comparer class for HashSet<T> objects.
/// </summary>
/// <typeparam name="T">The type of elements in the HashSet.</typeparam>
public class SetEquality<T> : EqualityComparer<HashSet<T>> where T : IEquatable<T>
{
    /// <summary>
    /// Determines whether two HashSets are equal. 
    /// </summary>
    /// <param name="x">The first HashSet to compare.</param>
    /// <param name="y">The second HashSet to compare.</param>
    /// <returns><c>true</c>, if both sets are equal; otherwise, <c>false</c>. </returns>
    public override bool Equals(HashSet<T>? x, HashSet<T>? y)
    {
        return y is not null && x is not null && x.SetEquals(y);
    }

    /// <summary> 
    /// Returns a hash code for the specified object. 
    /// </summary> 
    /// <param name="obj">The object for which a hash code is to be returned.</param> 
    /// <returns></returns> 		
    public override int GetHashCode(HashSet<T> obj) => (obj.Count(), typeof(T).GetHashCode()).GetHashCode();
}

/// <summary>
/// Represents an EqualityComparer for IEnumerable of type T, where T is IEquatable and IComparable.
/// </summary>
/// <typeparam name="T">The type of the elements of the IEnumerable.</typeparam>
public class SequenceEquality<T> : EqualityComparer<IEnumerable<T>> where T : IEquatable<T>, IComparable<T>
{
    /// <summary>
    /// Compares two IEnumerable objects for equality.
    /// </summary>
    /// <param name="x">The first IEnumerable object to compare.</param>
    /// <param name="y">The second IEnumerable object to compare.</param>
    /// <returns>True if the two objects are equal, false otherwise.</returns>
    public override bool Equals(IEnumerable<T>? x, IEnumerable<T>? y)
    {
        return y is not null && x is not null && x.SequenceEqual(y);
    }

    /// <summary>
    /// Gets the hash code for the given <see cref="IEnumerable{T}"/>.
    /// </summary>
    /// <param name="obj">The <see cref="IEnumerable{T}"/> to get the hash code for.</param>
    /// <returns>The hash code for the given <see cref="IEnumerable{T}"/>.</returns>
    public override int GetHashCode(IEnumerable<T> obj) => (obj.Count(), typeof(T).GetHashCode()).GetHashCode();
}
