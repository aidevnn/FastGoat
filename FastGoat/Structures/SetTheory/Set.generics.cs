namespace FastGoat.Structures.SetTheory;

public interface IElt<U> : IElt, IEquatable<U>, IComparable<U> where U : struct, IElt<U>
{
    IFSet<U> FSet { get; }
}

public interface IFSet<U> : IFSet, IEquatable<IFSet<U>> where U : struct, IElt<U>
{
    void AddElement(U e);
    bool Contains(U e);
    IEnumerable<U> AllElements();
}

public interface ISubSet<U> : IFSet<U> where U : struct, IElt<U>
{
    IFSet<U> UpperSet { get; }
    int CompareElt(U a, U b);
    bool SetEquals(ISubSet<U> set);
}

public class EltEquality<U> : EqualityComparer<U> where U : struct, IElt<U>
{
    public override bool Equals(U x, U y) => x.FSet.Equals(y.FSet) && x.HashCode == y.HashCode;
    public override int GetHashCode(U obj) => obj.HashCode;
}
