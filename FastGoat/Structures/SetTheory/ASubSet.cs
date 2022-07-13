namespace FastGoat.Structures.SetTheory;

public abstract class SubSet<U> : ISubSet<U> where U : struct, IElt<U>
{
    protected SubSet(IFSet<U> fSet)
    {
        UpperSet = fSet;
        Elts = new HashSet<U>(new EltEquality<U>());
        Infos = new DisplaySet<U>(this);
    }

    public IFSet<U> UpperSet { get; }

    protected HashSet<U> Elts { get; }

    public virtual int CompareElt(U a, U b) => a.CompareTo(b);
    public virtual void AddElement(U e)
    {
        UpperSet.AddElement(e);
        Elts.Add(e);
    }

    public bool Contains(U e) => Elts.Contains(e);
    public int Count => Elts.Count;
    public IEnumerable<U> AllElements() => Elts;

    public bool Equals(IFSet<U>? other) => other?.GetHashCode() == GetHashCode();
    public bool SetEquals(ISubSet<U> set)
    {
        return UpperSet.Equals(set.UpperSet) && Elts.SetEquals(set.AllElements());
    }

    public DisplaySet<U> Infos { get; protected set; }
    public void SetName(string name) => Infos.Name = name;
    public void DisplayHead() => Infos.DisplayHead();
    public void DisplayElements(string? name = null, string? infos = null)
    {
        Infos.SetDetails(name: name, infos: infos);
        Infos.DisplayElements();
    }
    public void DisplayTable(char symb, Func<U, U, U> fct) => Infos.DisplayTable(symb, fct);
}

public class EqSubSet<U> : EqualityComparer<SubSet<U>> where U : struct, IElt<U>
{
    public override bool Equals(SubSet<U>? x, SubSet<U>? y) => x is null || y is null ? false : x.SetEquals(y);
    public override int GetHashCode(SubSet<U> obj) => obj.UpperSet.GetHashCode();
}