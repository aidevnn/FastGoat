using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Perms;

public enum DisplayPerm
{
    Table,
    Cycles,
    TableComma,
    CyclesComma
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
        return Sn.N < 11
            ? Hash == other.Hash
            : Hash == other.Hash && Table.SequenceEqual(other.Table);
    }

    public override int GetHashCode()
    {
        return Hash;
    }

    public int[][] Orbits => IntExt.PermutationToCycles(Sn.N, Table);
    public int[] PermType => Orbits.Select(l => l.Length).Order().ToArray();
    public int Sgn => (-1).Pow(Orbits.Length);
    public string PermTypeStr => $"({PermType.Glue(" ")})";
    public string SgnStr => Sgn == 1 ? "(+)" : "(-)";

    public T[] Apply<T>(T[] ts)
    {
        if (ts.Length != Sn.N)
            throw new GroupException(GroupExceptionType.GroupDef);

        return Table.Select(i => ts[i]).ToArray();
    }

    public override string ToString()
    {
        if (Style == DisplayPerm.Table)
            return $"[{Table.Select(a => a + 1).Glue(" ")}]";
        if (Style == DisplayPerm.TableComma)
            return $"[{Table.Select(a => a + 1).Glue(", ")}]";

        var cycles = Orbits.Where(a => a.Length > 1).ToArray();
        if (Style == DisplayPerm.Cycles)
        {
            var strCycles = cycles.Select(a => $"({a.Select(b => b + 1).Glue(" ")})").Glue();
            return $"[{strCycles}]";
        }
        else
        {
            var strCycles = cycles.Select(a => $"({a.Select(b => b + 1).Glue(", ")})").Glue(", ");
            return $"[{strCycles}]";
        }
    }

    public override bool Equals(object? obj)
    {
        return obj is Perm perm && Equals(perm);
    }

    public static Perm operator *(Perm a, Perm b) => a.BaseGroup.Op(a, b);
    public static Perm operator ^(Perm a, int p) => a.BaseGroup.Times(a, p);
    public static bool operator ==(Perm a, Perm b) => a.Equals(b);
    public static bool operator !=(Perm a, Perm b) => !a.Equals(b);

    public static int CompareOrbits(Perm a, Perm b)
    {
        if (!a.Sn.Equals(b.Sn))
            throw new();

        var ca = a.Orbits.Where(e => e.Length > 1).ToArray();
        var cb = b.Orbits.Where(e => e.Length > 1).ToArray();
        ;

        var compNb = ca.Length.CompareTo(cb.Length);
        if (compNb != 0)
            return compNb;

        ca = ca.OrderBy(e => e, Comparer<int[]>.Create((e, f) => e.SequenceCompareTo(f))).ThenBy(e => e.Length).ToArray();
        cb = cb.OrderBy(e => e, Comparer<int[]>.Create((e, f) => e.SequenceCompareTo(f))).ThenBy(e => e.Length).ToArray();
        foreach (var (e, f) in ca.Zip(cb))
        {
            var compL = e.Length.CompareTo(f.Length);
            if (compL != 0)
                return compL;

            var compEF = e.SequenceCompareTo(f);
            if (compEF != 0)
                return compEF;
        }

        return 0;
    }

    public static Comparer<Perm> OrbitsComparer => Comparer<Perm>.Create(CompareOrbits);
}