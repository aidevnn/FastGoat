using System.Collections;
using FastGoat.Commons;
using FastGoat.UserGroup.Integers;
using Mono.Options;

namespace FastGoat.Structures.GenericGroup;

public readonly struct TableElt : IElt<TableElt>
{
    public int E { get; }
    public int GroupTableHash { get; }
    public string Str { get; }

    public TableElt(int groupTableHash, string str, int e)
    {
        (GroupTableHash, Str, E) = (groupTableHash, str, e);
        Hash = (GroupTableHash, E).GetHashCode();
    }

    public bool Equals(TableElt other) => GroupTableHash == other.GroupTableHash && E.Equals(other.E);

    public int CompareTo(TableElt other) => E.CompareTo(other.E);

    public int Hash { get; }

    public override int GetHashCode() => Hash;
    public override string ToString() => Str;
}

public readonly struct GroupTableBase : IGroup<TableElt>
{
    public Dictionary<(TableElt, TableElt), TableElt> OpTable { get; }
    public Dictionary<TableElt, TableElt> InvTable { get; }
    public int Order { get; }
    public TableElt NeutralElt { get; }

    public GroupTableBase(string name, int hash, TableElt neutral, Dictionary<(TableElt, TableElt), TableElt> opTable,
        Dictionary<TableElt, TableElt> invTable, TableElt[] gens)
    {
        (Name, Hash) = (name, hash);
        (NeutralElt, OpTable, InvTable) = (neutral, opTable, invTable);
        Elements = invTable.Keys.ToHashSet();
        PseudoGeneratos = gens;
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

    public TableElt Invert(TableElt e) => InvTable[e];

    public TableElt Op(TableElt e1, TableElt e2) => OpTable[(e1, e2)];

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{Name}(table)";

    public static (Dictionary<T, TableElt> idxTable, GroupTableBase gb) Create<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var elts = g.OrderBy(e => g.ElementsOrders[e]).ToArray();
        var idxTable = elts.Select((e, i) => (e, i)).ToDictionary(e => e.e, e => new TableElt(g.Hash, $"{e.e}", e.i));
        var invTable = elts.Select(e => (e, g.Invert(e))).ToDictionary(e => idxTable[e.e], e => idxTable[e.Item2]);
        var opTable = elts.Grid2D(elts).Select(e => (e, g.Op(e.t1, e.t2)))
            .ToDictionary(e => (idxTable[e.e.t1], idxTable[e.e.t2]), e => idxTable[e.Item2]);
        var gens = g.GetGenerators().Select(e => idxTable[e]).ToArray();
        return (idxTable, new(g.Name, g.Hash, idxTable[g.Neutral()], opTable, invTable, gens));
    }
}

public class GroupTable : ConcreteGroup<TableElt>
{
    private GroupTable(GroupTableBase tableBase) : base(tableBase, tableBase.GetGenerators().ToArray())
    {
    }

    public static (Dictionary<T, TableElt> idxTable, GroupTable gt) Create<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var (tb, gb) = GroupTableBase.Create(g);
        return (tb, new(gb));
    }
}