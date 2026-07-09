using System.Collections;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace Craft.Craft;

public struct CosetGroup<T> : IGroup<Coset<T>> where T : struct, IElt<T>
{
    public ConcreteGroup<T> G { get; }
    public ConcreteGroup<T> H { get; }
    public CosetType CosetType { get; }
    public Dictionary<T, Coset<T>> Map { get; }

    public CosetGroup(ConcreteGroup<T> g, ConcreteGroup<T> h, CosetType cosetType = CosetType.Both)
    {
        (G, H, CosetType) = (g, h, cosetType);
        Map = Group.Cosets(G, H, CosetType);
        var name = $"{G.NameParenthesis()}/{H.NameParenthesis()}";
        Name = CosetType switch
        {
            CosetType.Left => $"{name}(l)",
            CosetType.Right => $"{name}(r)",
            _ => name
        };

        Elements = Map.Values.ToHashSet();
        Hash = (G.Hash, H.Hash, cosetType).GetHashCode();
    }
    
    private HashSet<Coset<T>> Elements { get; }

    public IEnumerator<Coset<T>> GetEnumerator() => Elements.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<Coset<T>>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public Coset<T> this[params ValueType[] us] => throw new NotImplementedException();

    public IEnumerable<Coset<T>> GetElements() => Elements;

    public IEnumerable<Coset<T>> GetGenerators()
    {
        yield return Neutral();
    }

    public Coset<T> Neutral() => new(G, H, G.Neutral(), CosetType);

    public Coset<T> Invert(Coset<T> e) => Map[G.Invert(e.X)];

    public Coset<T> Op(Coset<T> e1, Coset<T> e2) => CosetType == CosetType.Left ? Map[G.Op(e1.X, e2.X)] : Map[G.Op(e2.X, e1.X)];

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}