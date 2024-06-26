using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;

namespace FastGoat.UserGroup.Perms;

public struct Sn : IGroup<Perm>
{
    public int Hash { get; }
    public string Name { get; }
    public int N { get; }
    readonly int[] _cache;

    public Sn()
    {
        N = Hash = 2;
        _cache = new int[2];
        Name = $"S{N}";
    }

    public Sn(int n)
    {
        if (n < 1)
            throw new GroupException(GroupExceptionType.GroupDef);

        N = Hash = n;
        _cache = new int[n];
        Name = $"S{N}";
    }

    public IEnumerable<Perm> GetGenerators()
    {
        if (N > 1)
        {
            yield return this[(1, 2)];
            yield return ComposesCycles(new Tuple2Array(Enumerable.Range(1, N).ToArray()));
        }
        else
            yield return this[1];
    }

    public IEnumerable<Perm> GetElements()
    {
        return IntExt.GetPermutations(N).Select(CreateElement);
    }

    public Perm Neutral() => new Perm(this);

    public Perm Invert(Perm e)
    {
        if (!Equals(e.BaseGroup))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var hash = IntExt.InvertPermutation(e.Table, _cache);
        return new Perm(this, _cache, hash);
    }

    public Perm Op(Perm e1, Perm e2)
    {
        if (!Equals(e1.BaseGroup) || !Equals(e2.BaseGroup))
            throw new GroupException(GroupExceptionType.BaseGroup);

        var hash = IntExt.ComposePermutation(e1.Table, e2.Table, _cache);
        return new Perm(this, _cache, hash);
    }

    public Perm CreateElement(params int[] table)
    {
        var t0 = table.Select(i => i - 1).ToArray();
        if (!IntExt.CheckTable(N, t0))
            throw new GroupException(GroupExceptionType.GroupDef);

        var hash = IntExt.GenHash(N, t0);
        var p = new Perm(this, t0, hash);
        return p;
    }

    public Perm CreateElementTable(params int[] table)
    {
        if (!IntExt.CheckTable(N, table))
            throw new GroupException(GroupExceptionType.GroupDef);

        var hash = IntExt.GenHash(N, table);
        var p = new Perm(this, table, hash);
        return p;
    }

    public Perm Cycle(params int[] cycle) => ComposesCycles(new Tuple2Array(cycle));

    public Perm ComposesCycles(params Tuple2Array[] cycles)
    {
        Neutral().Table.CopyTo(_cache, 0);
        foreach (var e in cycles)
        {
            var cycle = e.Table.Select(i => i - 1).ToArray();
            if (!IntExt.CheckCycle(N, cycle))
                throw new GroupException(GroupExceptionType.GroupDef);

            IntExt.ApplyCycle(_cache, cycle);
        }

        var hash = IntExt.GenHash(N, _cache);
        return new Perm(this, _cache, hash);
    }

    public bool Equals(IGroup<Perm>? other) => other?.Hash == Hash;

    public Perm this[params ValueType[] us]
    {
        get
        {
            var cycles = us.Select(u => (Tuple2Array)u).ToArray();
            if (cycles.Any(c => c.Table.Length == 0))
                throw new GroupException(GroupExceptionType.GroupDef);

            if (cycles.All(c => c.Table.Length == 1) && cycles.Length == N)
            {
                var table = cycles.SelectMany(c => c.Table).ToArray();
                return CreateElement(table);
            }

            if (cycles.Any(c => c.Table.Length == 1))
                throw new GroupException(GroupExceptionType.GroupDef);

            return ComposesCycles(cycles);
        }
    }

    public IEnumerator<Perm> GetEnumerator() => GetElements().GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => GetElements().GetEnumerator();
    public override string ToString() => Name;
}
