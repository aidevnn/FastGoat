using System.Runtime.CompilerServices;

namespace FastGoat.Commons;

public class Tuple2Array
{
    public int[] Table { get; private set; }

    public Tuple2Array(params int[] table)
    {
        Table = table.ToArray();
    }

    public override string ToString() => $"[{Table.Glue(" ")}]";

    public static implicit operator Tuple2Array(ValueType v)
    {
        if (v is int v1)
            return new(v1);

        if (v is ITuple vx && Enumerable.Range(0, vx.Length).All(i => vx[i] is int))
            return new(Enumerable.Range(0, vx.Length).Select(i => (int)(vx[i] ?? 0)).ToArray());

        return new();
    }

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