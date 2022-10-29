using FastGoat.Commons;
using FastGoat.Theory;

namespace FastGoat.UserGroup.Perms;

public enum DisplayPerm
{
    Table,
    Cycles
}

public struct Perm : IElt<Perm>
{
    public static DisplayPerm Style { get; set; } = DisplayPerm.Cycles;
    public int Hash { get; }
    public int[] Table { get; }
    public IGroup<Perm> BaseGroup { get; }
    public Sn Sn { get; }

    public Perm()
    {
        BaseGroup = Sn = new Sn(2);
        Table = new[] { 0, 1 };
        Hash = IntExt.GenHash(2, Table);
    }

    public Perm(Sn sn)
    {
        BaseGroup = Sn = sn;
        Table = sn.N.Range();
        Hash = IntExt.GenHash(sn.N, Table);
    }

    public Perm(Sn sn, int[] arr, int hash)
    {
        BaseGroup = Sn = sn;
        Table = arr.ToArray();
        Hash = hash;
    }

    public int CompareTo(Perm other)
    {
        if (!BaseGroup.Equals(other.BaseGroup))
            throw new GroupException(GroupExceptionType.BaseGroup);

        return Table.SequenceCompareTo(other.Table);
    }

    public bool Equals(Perm other)
    {
        return Hash == other.Hash;
    }

    public override int GetHashCode()
    {
        return Hash;
    }

    public override string ToString()
    {
        if (Style == DisplayPerm.Table)
            return $"[{Table.Select(a => a + 1).Glue(" ")}]";

        var orbits = IntExt.PermutationToCycles(Sn.N, Table);
        var cycles = orbits.Where(a => a.Length > 1).ToArray();
        var strCycles = cycles.Select(a => $"({a.Select(b => b + 1).Glue(" ")})").Glue();
        return $"[{strCycles}]";
    }

    public override bool Equals(object? obj)
    {
        return obj is Perm perm && Equals(perm);
    }

    public static Perm operator *(Perm a, Perm b) => a.BaseGroup.Op(a, b);
    public static Perm operator ^(Perm a, int p) => a.BaseGroup.Times(a, p);
    public static bool operator ==(Perm a, Perm b) => a.Equals(b);
    public static bool operator !=(Perm a, Perm b) => !a.Equals(b);
}