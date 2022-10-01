namespace FastGoat;

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

        if (v is ValueTuple<int, int> v2)
        {
            (int e1, int e2) p = v2;
            return new(p.e1, p.e2);
        }

        if (v is ValueTuple<int, int, int> v3)
        {
            (int e1, int e2, int e3) p = v3;
            return new(p.e1, p.e2, p.e3);
        }

        if (v is ValueTuple<int, int, int, int> v4)
        {
            (int e1, int e2, int e3, int e4) p = v4;
            return new(p.e1, p.e2, p.e3, p.e4);
        }

        if (v is ValueTuple<int, int, int, int, int> v5)
        {
            (int e1, int e2, int e3, int e4, int e5) p = v5;
            return new(p.e1, p.e2, p.e3, p.e4, p.e5);
        }

        if (v is ValueTuple<int, int, int, int, int, int> v6)
        {
            (int e1, int e2, int e3, int e4, int e5, int e6) p = v6;
            return new(p.e1, p.e2, p.e3, p.e4, p.e5, p.e6);
        }

        if (v is ValueTuple<int, int, int, int, int, int, int> v7)
        {
            (int e1, int e2, int e3, int e4, int e5, int e6, int e7) p = v7;
            return new(p.e1, p.e2, p.e3, p.e4, p.e5, p.e6, p.e7);
        }

        return new();
    }
}