using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public readonly struct Coset<T> : ICoset<T, Coset<T>> where T : struct, IElt<T>
{
    public Coset(ConcreteGroup<T> g, ConcreteGroup<T> h, CosetType cosetType = CosetType.Both)
    {
        CosetType = cosetType;
        G = g;
        H = h;
        X = H.Neutral();
        Hash = (g.Hash, h.Hash, X.Hash).GetHashCode();
    }

    public Coset(ConcreteGroup<T> g, ConcreteGroup<T> h, T x, CosetType cosetType = CosetType.Both)
    {
        CosetType = cosetType;
        G = g;
        H = h;
        X = x;
        Hash = (g.Hash, h.Hash, X.Hash).GetHashCode();
    }
    
    public CosetType CosetType { get; }

    public T X { get; }
    public ConcreteGroup<T> G { get; }
    public ConcreteGroup<T> H { get; }
    public bool Equals(Coset<T> other) => X.Equals(other.X);

    public int CompareTo(Coset<T> other) => X.CompareTo(other.X);

    private IEnumerable<T> Elements
    {
        get
        {
            List<T> set = new();
            foreach (var h in H)
            {
                if (CosetType != CosetType.Right)
                    set.Add(G.Op(X, h));
                else
                    set.Add(G.Op(h, X));
            }

            return set.Ascending();
        }
    }

    public IEnumerator<T> GetEnumerator() => Elements.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var hName = H.Name;
        if (CosetType != CosetType.Right)
            return $"{X}({hName})";
        else
            return $"({hName}){X}";
    }
}
