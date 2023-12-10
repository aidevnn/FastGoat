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
    public Dictionary<(TableElt, TableElt), TableElt> OpTable { get; }
    public Dictionary<TableElt, TableElt> InvTable { get; }
    public int Order { get; }
    public TableElt NeutralElt { get; }
    private dynamic BaseGroup { get; }
    private Dictionary<int, TableElt> IdxTable { get; }

    public GroupTableBase(dynamic group, Dictionary<int, TableElt> idxTable, TableElt[] gens)
    {
        (Name, Hash) = (group.Name, (group.Hash, "table").GetHashCode());
        BaseGroup = group;
        IdxTable = idxTable;
        PseudoGeneratos = gens;
        var og = IdxTable.Count;
        OpTable = new(og * og);
        InvTable = new(og);
        Elements = IdxTable.Values.ToHashSet();
        NeutralElt = IdxTable[group.Neutral().GetHashCode()];
        Order = Elements.Count;
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
        if (!InvTable.ContainsKey(e))
        {
            var e0 = BaseGroup.Invert(e.E).GetHashCode();
            var ei = InvTable[e] = IdxTable[e0];
            return ei;
        }
        
        return InvTable[e];
    }

    public TableElt Op(TableElt e1, TableElt e2)
    {
        var e1e2 = (e1, e2);
        if (!OpTable.ContainsKey(e1e2))
        {
            var e12Hash = BaseGroup.Op(e1.E, e2.E).GetHashCode();
            var e12 = OpTable[e1e2] = IdxTable[e12Hash];
            return e12;
        }

        return OpTable[e1e2];
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
    private GroupTable(GroupTableBase tableBase) : base(tableBase, tableBase.GetGenerators().ToArray())
    {
    }

    public static GroupTable Create<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return new(GroupTableBase.Create(g));
    }
}