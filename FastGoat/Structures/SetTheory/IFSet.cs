namespace FastGoat.Structures.SetTheory;

public interface IFSet<U> : IEquatable<IFSet<U>> where U : struct, IElt<U>
{
    void AddElement(U e);
    bool Contains(U e);
    IEnumerable<U> AllElements { get; }
    int Count { get; }
}

public interface ISubSet<U> : IFSet<U> where U : struct, IElt<U>
{
    IFSet<U> UpperSet { get; }
    int CompareElt(U a, U b);
    bool SetEquals(ISubSet<U> set);
}
