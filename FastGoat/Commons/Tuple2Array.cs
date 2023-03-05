using System.Runtime.CompilerServices;

namespace FastGoat.Commons;

/// <summary>
/// Extension class for implicit conversion of integers tuples to integers array. 
/// </summary>
public class Tuple2Array
{
    /// <summary>
    /// Gets or sets the Table property.
    /// </summary>
    public int[] Table { get; private set; }

    /// <summary>
    /// Creates a new instance of Tuple2Array with the given parameters.
    /// </summary>
    /// <param name="table">An array of integers representing the elements in the tuple.</param>
    public Tuple2Array(params int[] table)
    {
        Table = table.ToArray();
    }

    /// <summary>
    /// Returns a string representation of the Table object.
    /// </summary>
    /// <returns>A string containing the elements of the Table object, separated by spaces.</returns>
    public override string ToString() => $"[{Table.Glue(" ")}]";

    /// <summary>
    /// Implicitly converts a ValueType to a Tuple2Array.
    /// </summary>
    /// <remarks>The ValueType can be an integer or a tuple of integers</remarks>
    /// <param name="v">The ValueType to convert.</param>
    /// <returns>A new Tuple2Array.</returns>
    public static implicit operator Tuple2Array(ValueType v)
    {
        if (v is int v1)
            return new(v1);

        if (v is ITuple vx && Enumerable.Range(0, vx.Length).All(i => vx[i] is int))
            return new(Enumerable.Range(0, vx.Length).Select(i => (int)(vx[i] ?? 0)).ToArray());

        return new();
    }

    /// <summary>
    /// Converts a ValueType to a Tuple2Array.
    /// </summary>
    /// <remarks>The ValueType can be an integer or a tuple of integers or a tuple of tuples of integer</remarks>
    /// <param name="v">The ValueType to convert.</param>
    /// <returns>A new array of Tuple2Array.</returns>
    public static Tuple2Array[] ComplexTuples(ValueType v)
    {
        var c0 = (Tuple2Array)v;
        if (c0.Table.Length > 0)
            return new[] { c0 };

        if (v is ITuple t)
        {
            var tuples = Enumerable.Range(0, t.Length)
                .Select(i => t[i] as ValueType ?? new Tuple2Array())
                .Where(e => e.Table.Length > 0)
                .ToArray();

            return tuples;
        }

        return Array.Empty<Tuple2Array>();
    }
}