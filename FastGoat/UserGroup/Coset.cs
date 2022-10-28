using System.Collections;
using FastGoat;

namespace FastGoat.UserGroup;

public readonly struct Coset<T> : ILeftCoset<T, Coset<T>> where T : struct, IElt<T>
{
    public Coset(Quotient<T> lQuo)
    {
        Quotient = lQuo;
        X = lQuo.H.Neutral();
        Hash = (lQuo.Hash, X.Hash).GetHashCode();
    }

    public Coset(Quotient<T> lQuo, T x)
    {
        Quotient = lQuo;
        X = x;
        Hash = (lQuo.Hash, X.Hash).GetHashCode();
    }

    public T X { get; }
    public Quotient<T> Quotient { get; }
    public bool Equals(Coset<T> other) => Hash == other.Hash;

    public int CompareTo(Coset<T> other) => X.CompareTo(other.X);

    public IEnumerable<T> xH
    {
        get
        {
            List<T> set = new();
            foreach (var h in Quotient.H)
            {
                set.Add(Quotient.G.Op(X, h));
            }

            return set.Ascending();
        }
    }

    public int Hash { get; }
    public IGroup<Coset<T>> BaseGroup => Quotient;
    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var hName = Quotient.H.Name;
        return $"{X}({hName})";
    }
}