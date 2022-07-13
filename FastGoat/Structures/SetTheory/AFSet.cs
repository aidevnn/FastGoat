namespace FastGoat.Structures.SetTheory;

public abstract class FSet<U> : IFSet<U> where U : struct, IElt<U>
{
    protected FSet()
    {
        Elts = new HashSet<U>(new EltEquality<U>());
    }

    protected HashSet<U> Elts { get; }

    public void AddElement(U e) => Elts.Add(e);
    public bool Contains(U e) => Elts.Contains(e);
    public int Count => Elts.Count;
    public IEnumerable<U> AllElements() => Elts;

    public bool Equals(IFSet<U>? other) => other?.GetHashCode() == GetHashCode();
}