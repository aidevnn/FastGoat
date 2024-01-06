using System.Collections;
using FastGoat.Commons;
using FastGoat.UserGroup.Integers;
using Mono.Options;

namespace FastGoat.Structures.GenericGroup;

public readonly struct TableElt : IElt<TableElt>
{
    public dynamic E { get; }

    public TableElt(dynamic e)
    {
        E = e;
        Hash = E.GetHashCode();
    }

    public bool Equals(TableElt other) => Hash == other.Hash;

    public int CompareTo(TableElt other) => E.CompareTo(other.E);

    public int Hash { get; }

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{E}";
}

public readonly struct GroupTableBase : IGroup<TableElt>
{
    public TableElt NeutralElt { get; }
    private dynamic OriginGroup { get; }
    private Dictionary<int, TableElt> IdxTable { get; }

    public GroupTableBase(dynamic group, Dictionary<int, TableElt> idxTable, TableElt[] gens)
    {
        (Name, Hash) = (group.Name, (group.Hash, "table").GetHashCode());
        OriginGroup = group;
        IdxTable = idxTable;
        PseudoGeneratos = gens;
        Elements = IdxTable.Values.ToHashSet();
        NeutralElt = IdxTable[group.Neutral().GetHashCode()];
    }

    public HashSet<TableElt> Elements { get; }
    public TableElt[] PseudoGeneratos { get; }

    public IEnumerator<TableElt> GetEnumerator() => Elements.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<TableElt>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public TableElt this[params ValueType[] us] => throw new NotImplementedException();

    public IEnumerable<TableElt> GetElements() => Elements;

    public IEnumerable<TableElt> GetGenerators() => PseudoGeneratos;

    public TableElt Neutral() => NeutralElt;

    public TableElt Invert(TableElt e)
    {
        var e0 = OriginGroup.Invert(e.E).GetHashCode();
        return IdxTable[e0];
    }

    public TableElt Op(TableElt e1, TableElt e2)
    {
        var e12Hash = OriginGroup.Op(e1.E, e2.E).GetHashCode();
        return IdxTable[e12Hash];
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{Name.WithParenthesis()}(table)";

    public static GroupTableBase Create<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var hashes = g.Select(e => e.Hash).ToHashSet();
        if (hashes.Count != g.Count())
            throw new("Hash collisions");

        var idxTable = g.Order().ToDictionary(e => e.Hash, e => new TableElt(e));
        var gens = g.GetGenerators().Select(e => new TableElt(e)).ToArray();
        return new GroupTableBase(g, idxTable, gens);
    }
}

public class GroupTable : ConcreteGroup<TableElt>
{
    private GroupTableBase TableBase { get; }
    private GroupTable(GroupTableBase tableBase) : base(tableBase, tableBase.GetGenerators().ToArray())
    {
        TableBase = tableBase;
    }

    public static GroupTable Create<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return new(GroupTableBase.Create(g));
    }
}
