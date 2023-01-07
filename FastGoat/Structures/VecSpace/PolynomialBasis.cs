using FastGoat.Commons;

namespace FastGoat.Structures.VecSpace;

public readonly struct PolynomialBasis<K, T> : IElt<PolynomialBasis<K, T>>
    where K : struct, IElt<K>, IRingElt<K>, IFieldElt<K>
    where T : struct, IElt<T>
{
    public Polynomial<K, T>[] Basis { get; }
    public K[] LC { get; }
    public Monom<T>[] LM { get; }
    public Polynomial<K, T>[] LT { get; }
    public Indeterminates<T> Indeterminates { get; }

    public PolynomialBasis()
    {
        throw new ArgumentException();
    }

    public PolynomialBasis(Indeterminates<T> indeterminates, params Polynomial<K, T>[] basis)
    {
        Indeterminates = indeterminates;
        Basis = basis;
        var n = basis.Length;
        LC = new K[n];
        LM = new Monom<T>[n];
        LT = new Polynomial<K, T>[n];
        Hash = 0;
        for (int i = 0; i < n; ++i)
        {
            var f = Basis[i];
            if (!f.Indeterminates.Equals(Indeterminates))
                throw new ArgumentException();
            
            (LC[i], LM[i], LT[i]) = f.LeadingDetails;
            Hash = (Hash, f.Hash).GetHashCode();
        }
    }

    public int Hash { get; }
    public Polynomial<K, T> Rem(Polynomial<K, T> f)
    {
        var f0 = f;
        foreach (var p in Basis)
            f0 = f0.Div(p).rem;

        return f0;
    }

    public bool Equals(PolynomialBasis<K, T> other) => Basis.SequenceEqual(other.Basis);

    public int CompareTo(PolynomialBasis<K, T> other) => Basis.SequenceCompareTo(other.Basis);

    public override int GetHashCode() => Hash;
    public override string ToString() => $"[{Basis.Glue("; ")}]";
}