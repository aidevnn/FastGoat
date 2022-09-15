namespace FastGoat;

public struct Sn : IGroup<Perm>
{
    public int Hash { get; }
    public int N { get; }
    readonly int[] cache;
    public Sn()
    {
        throw new Exception();
    }
    public Sn(int n)
    {
        // if (n < 2 || n > 8)
        //     throw new Exception("Seventh Sky is the limit.");

        N = Hash = n;
        cache = new int[n];
    }
    public Perm Neutral() => new Perm(this);
    public Perm Invert(Perm a)
    {
        var hash = IntExt.InvertPermutation(a.Table, cache);
        return new Perm(this, cache, hash);
    }
    public Perm Op(Perm a, Perm b)
    {
        var hash = IntExt.ComposePermutation(a.Table, b.Table, cache);
        return new Perm(this, cache, hash);
    }
    public Perm CreateElement(params int[] vs)
    {
        vs.Add(-1);
        if (!IntExt.CheckTable(N, vs))
            return Neutral();

        var hash = IntExt.GenHash(N, vs);
        var p = new Perm(this, vs, hash);
        return p;
    }
    public Perm Cycle(params int[] vs)
    {
        vs.Add(-1);
        if (!IntExt.CheckCycle(N, vs))
            return Neutral();

        Neutral().Table.CopyTo(cache, 0);
        IntExt.ApplyCycle(cache, vs);
        var hash = IntExt.GenHash(N, cache);
        return new Perm(this, cache, hash);
    }

    public Perm KCycle(int count) => Cycle(count.Range(1));
    public Perm KCycle(int start, int count) => Cycle(count.Range(start));

    public Perm Cycle(Tuple2Array cycle) => Cycle(cycle.Table);
    public Perm ComposesCycles(params Tuple2Array[] cycles)
    {
        Neutral().Table.CopyTo(cache, 0);
        foreach (var e in cycles)
        {
            e.Table.Add(-1);
            if (!IntExt.CheckCycle(N, e.Table)) return Neutral();

            IntExt.ApplyCycle(cache, e.Table);
        }

        var hash = IntExt.GenHash(N, cache);
        return new Perm(this, cache, hash);
    }

    public Perm C(params int[] vs) => Cycle(vs);
    public Perm C(params Tuple2Array[] cycles) => ComposesCycles(cycles);

    public bool Equals(IGroup<Perm>? other) => other?.Hash == Hash;

    public Perm this[int k]
    {
        get
        {
            var arr = IntExt.GetPermutation(N, k).ToArray();
            return CreateElement(arr);
        }
    }

    public override string ToString() => $"S{N}";
}
