using FastGoat;

namespace FastGoat.UserGroup;

public readonly struct Representative<T> : IElt<Representative<T>> where T : struct, IElt<T>
{
    public Representative(Quotient<T> lQuo)
    {
        Quotient = lQuo;
        X = lQuo.H.Neutral();
        Hash = (lQuo.Hash, X.Hash).GetHashCode();
    }

    public Representative(Quotient<T> lQuo, T x)
    {
        Quotient = lQuo;
        X = x;
        Hash = (lQuo.Hash, X.Hash).GetHashCode();
    }

    public T X { get; }
    public Quotient<T> Quotient { get; }
    public bool Equals(Representative<T> other) => Hash == other.Hash;

    public int CompareTo(Representative<T> other) => X.CompareTo(other.X);

    public int Hash { get; }
    public IGroup<Representative<T>> BaseGroup => Quotient;
    public override int GetHashCode() => Hash;
    public override string ToString()
    {
        var hName = Quotient.H.Name;
        return hName.Contains(' ') ? $"{X}_({hName})" : $"{X}_{hName}";
    }
}