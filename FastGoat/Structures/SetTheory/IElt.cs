namespace FastGoat.Structures.SetTheory;

public interface IElt<U> : IEquatable<U>, IComparable<U> where U : struct, IElt<U>
{
    IFSet<U> FSet { get; }
    int HashCode { get; }
    int GetHashCode();
    string ToString();
}

public class EltEquality<U> : EqualityComparer<U> where U : struct, IElt<U>
{
    public override bool Equals(U x, U y) => x.FSet.Equals(y.FSet) && x.Equals(y);
    public override int GetHashCode(U obj) => obj.GetHashCode();
}
