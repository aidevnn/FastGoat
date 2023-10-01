using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures.CartesianProduct;

namespace FastGoat.Structures.GenericGroup;

public struct MapGroupBase<T1, T2> : IGroup<MapElt<T1, T2>> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    public MapGroupBase(ConcreteGroup<T1> G1, ConcreteGroup<T2> G2)
    {
        (this.G2, Domain) = (G2, G1);
        Hash = (G1.Count(), G2.Count()).GetHashCode();
        Name = $"MapGroups({Domain},{this.G2})";
    }
    
    public ConcreteGroup<T1> Domain { get; }
    public ConcreteGroup<T2> G2 { get; }
    public IEnumerator<MapElt<T1, T2>> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<MapElt<Ep<T1>, T2>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public MapElt<T1, T2> this[params ValueType[] us] => throw new NotImplementedException();

    public IEnumerable<MapElt<T1, T2>> GetElements()
    {
        yield return Neutral();
    }

    public IEnumerable<MapElt<T1, T2>> GetGenerators() 
    {
        yield return Neutral();
    }

    public MapElt<T1, T2> Neutral()
    {
        var nN = G2.Neutral();
        var map = Domain.ToDictionary(g0 => g0, _ => nN);
        return new MapElt<T1, T2>(Domain, G2, map);
    }

    public MapElt<T1, T2> Invert(MapElt<T1, T2> e)
    {
        var map = new Dictionary<T1, T2>();
        foreach (var v in Domain)
            map[v] = G2.Invert(e[v]);

        return new MapElt<T1, T2>(Domain, G2, map);
    }

    public MapElt<T1, T2> Op(MapElt<T1, T2> e1, MapElt<T1, T2> e2)
    {
        var map = new Dictionary<T1, T2>();
        foreach (var v in Domain)
            map[v] = G2.Op(e1[v], e2[v]);

        return new MapElt<T1, T2>(Domain, G2, map);
    }

    public bool Equals(MapGroupBase<T2, T1> other) => other.Hash == Hash;

    public bool Equals(IGroup<MapElt<T1, T2>>? other) => other?.Hash == Hash;

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}