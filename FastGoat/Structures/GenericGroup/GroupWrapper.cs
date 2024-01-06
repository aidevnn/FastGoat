using System.Collections;
using FastGoat.Commons;
using FastGoat.UserGroup.Integers;
using Mono.Options;

namespace FastGoat.Structures.GenericGroup;

public readonly struct WElt : IElt<WElt>
{
    public dynamic E { get; }

    public WElt(dynamic e)
    {
        E = e;
        Hash = E.GetHashCode();
    }

    public bool Equals(WElt other) => Hash == other.Hash;

    public int CompareTo(WElt other) => E.CompareTo(other.E);

    public int Hash { get; }

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{E}";
}

public readonly struct GroupWrapperBase : IGroup<WElt>
{
    public WElt NeutralElt { get; }
    private dynamic OriginGroup { get; }
    private Dictionary<int, WElt> IdxTable { get; }

    public GroupWrapperBase(dynamic group, Dictionary<int, WElt> idxTable, WElt[] gens)
    {
        (Name, Hash) = (group.Name, (group.Hash, "table").GetHashCode());
        OriginGroup = group;
        IdxTable = idxTable;
        PseudoGeneratos = gens;
        Elements = IdxTable.Values.ToHashSet();
        NeutralElt = IdxTable[group.Neutral().GetHashCode()];
    }

    public HashSet<WElt> Elements { get; }
    public WElt[] PseudoGeneratos { get; }

    public IEnumerator<WElt> GetEnumerator() => Elements.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<WElt>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public WElt this[params ValueType[] us] => throw new NotImplementedException();

    public IEnumerable<WElt> GetElements() => Elements;

    public IEnumerable<WElt> GetGenerators() => PseudoGeneratos;

    public WElt Neutral() => NeutralElt;

    public WElt Invert(WElt e)
    {
        var e0 = OriginGroup.Invert(e.E).GetHashCode();
        return IdxTable[e0];
    }

    public WElt Op(WElt e1, WElt e2)
    {
        var e12Hash = OriginGroup.Op(e1.E, e2.E).GetHashCode();
        return IdxTable[e12Hash];
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => $"{Name.WithParenthesis()}(table)";

    public static GroupWrapperBase Create<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        var hashes = g.Select(e => e.Hash).ToHashSet();
        if (hashes.Count != g.Count())
            throw new("Hash collisions");

        var idxTable = g.Order().ToDictionary(e => e.Hash, e => new WElt(e));
        var gens = g.GetGenerators().Select(e => new WElt(e)).ToArray();
        return new GroupWrapperBase(g, idxTable, gens);
    }
}

public class GroupWrapper : ConcreteGroup<WElt>
{
    private GroupWrapperBase WrapperBase { get; }
    private GroupWrapper(GroupWrapperBase wrapperBase) : base(wrapperBase, wrapperBase.GetGenerators().ToArray())
    {
        WrapperBase = wrapperBase;
    }

    public static GroupWrapper Create<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        return new(GroupWrapperBase.Create(g));
    }
}
