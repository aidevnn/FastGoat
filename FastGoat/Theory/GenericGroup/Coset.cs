using System.Collections;
using FastGoat.Commons;

namespace FastGoat.Theory.GenericGroup;

public readonly struct Coset<T> : ILeftCoset<T, Coset<T>> where T : struct, IElt<T>
{
    public Coset(ConcreteGroup<T> g, ConcreteGroup<T> h)
    {
        G = g;
        H = h;
        X = H.Neutral();
        Hash = (g.Hash, h.Hash, X.Hash).GetHashCode();
    }

    public Coset(ConcreteGroup<T> g, ConcreteGroup<T> h, T x)
    {
        G = g;
        H = h;
        X = x;
        Hash = (g.Hash, h.Hash, X.Hash).GetHashCode();
    }

    public T X { get; }
    public ConcreteGroup<T> G { get; }
    public ConcreteGroup<T> H { get; }
    public bool Equals(Coset<T> other) => Hash == other.Hash;

    public int CompareTo(Coset<T> other) => X.CompareTo(other.X);

    private IEnumerable<T> xH
    {
        get
        {
            List<T> set = new();
            foreach (var h in H)
                set.Add(G.Op(X, h));

            return set.Ascending();
        }
    }

    public IEnumerator<T> GetEnumerator() => xH.GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public int Hash { get; }
    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var hName = H.Name;
        return $"{X}({hName})";
    }
}