using System.Collections;

namespace FastGoat;

public partial class ConcreteGroup<T> : IConcreteGroup<T> where T : struct, IElt<T>
{
    public int Hash { get; }
    public IGroup<T> BaseGroup { get; }
    public IConcreteGroup<T> ControlGroup { get; protected set; }
    protected HashSet<T> elements { get; set; }
    protected Dictionary<T, int> elementOrder { get; set; }
    public ConcreteGroup(IGroup<T> baseGroup)
    {
        BaseGroup = baseGroup;
        ControlGroup = this;
        Hash = this.GetHashCode();

        elements = new() { baseGroup.Neutral() };
        elementOrder = new() { [baseGroup.Neutral()] = 1 };
    }

    public IEnumerator<T> GetEnumerator() => elements.GetEnumerator();
    IEnumerator IEnumerable.GetEnumerator() => elements.GetEnumerator();
    public bool Equals(IConcreteGroup<T>? other) => other?.Hash == Hash;
    public GroupType groupType { get; protected set; }
    public int GetOrderOf(T e) => elementOrder[e];
    public virtual T Invert(T a) => BaseGroup.Invert(a);
    public virtual T Neutral() => BaseGroup.Neutral();
    public virtual T Op(T a, T b) => BaseGroup.Op(a, b);
    public T this[int k] => elements.ElementAt(k);
    public IEnumerable<T> SortByOrder(IEnumerable<T> ts) => ts.OrderBy(GetOrderOf).ThenAscending();
    public bool GroupEqual(IEnumerable<T> ts) => elements.SetEquals(ts);
    public bool IsSubGroupOf(IConcreteGroup<T> gr) => elements.IsSubsetOf(gr);
    public bool IsSuperGroupOf(IConcreteGroup<T> gr) => elements.IsSupersetOf(gr);
    public bool IsProperSubGroupOf(IConcreteGroup<T> gr) => elements.IsProperSubsetOf(gr);
    public bool IsProperSuperGroupOf(IConcreteGroup<T> gr) => elements.IsProperSupersetOf(gr);
    public void MinimalInfos()
    {
        Console.WriteLine($"Total:{this.Count()}");
        Console.WriteLine("Hashes => this:{0,12} SuperGroup:{1,12} BaseGroup:{2,12}", Hash, ControlGroup.Hash, BaseGroup.Hash);
        Console.WriteLine();
    }
}