namespace FastGoat;

public partial struct Perm : IElt<Perm>
{
    public static DisplayPerm Style = DisplayPerm.Cycles;
    public int Hash { get; }
    public int[] Table { get; }
    public IGroup<Perm> Group { get; }
    public Perm()
    {
        throw new Exception();
    }
    public Perm(Sn sn)
    {
        Group = sn;
        Table = sn.N.Range();
        Hash = IntExt.GenHash(sn.N, Table);
    }

    public Perm(Sn sn, int[] arr, int hash)
    {
        Group = sn;
        Table = arr.ToArray();
        Hash = hash;
    }

    public int CompareTo(Perm other) => Table.SequenceCompare(other.Table);
    public bool Equals(Perm other) => Hash == other.Hash;
    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        if (Style == DisplayPerm.Table)
            return $"[{Table.Select(a => a + 1).Glue(" ")}]";

        var orbits = Table.Orbits();
        var Cycles = orbits.Where(a => a.Count > 1).Select(a => a.ToArray()).ToArray();
        var strCycles = Cycles.Select(a => string.Format("({0})", a.Select(b => b + 1).Glue(" "))).Glue();
        return $"[{strCycles}]";
    }

    public static Perm operator *(Perm a, Perm b) => a.Group.Op(a, b);
    public static Perm operator ^(Perm a, int k)
    {
        var g = a.Group;
        if (k == 0)
            return g.Neutral();

        if (k < 0)
        {
            var ai = g.Invert(a);
            return Enumerable.Repeat(ai, -k).Aggregate((e0, e1) => g.Op(e0, e1));
        }

        return Enumerable.Repeat(a, k).Aggregate((e0, e1) => g.Op(e0, e1));
    }

    public static implicit operator Perm((Sn sn, Tuple2Array c0) p) => p.sn.ComposesCycles(p.c0);
    public static implicit operator Perm((Sn sn, Tuple2Array c0, Tuple2Array c1) p) => p.sn.ComposesCycles(p.c0, p.c1);
    public static implicit operator Perm((Sn sn, Tuple2Array c0, Tuple2Array c1, Tuple2Array c2) p) => p.sn.ComposesCycles(p.c0, p.c1, p.c2);
    public static implicit operator Perm((Sn sn, Tuple2Array c0, Tuple2Array c1, Tuple2Array c2, Tuple2Array c3) p) => p.sn.ComposesCycles(p.c0, p.c1, p.c2, p.c3);
}
