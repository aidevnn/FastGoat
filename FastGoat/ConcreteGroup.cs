using System.Collections;
using System.Collections.ObjectModel;

namespace FastGoat;

public class ConcreteGroup<T> : IConcreteGroup<T> where T : struct, IElt<T>
{
    public ConcreteGroup(string name, IGroup<T> g, bool singleton = false)
    {
        Name = name;
        Hash = new HashCode().ToHashCode();
        var ne = g.Neutral();
        BaseGroup = ne.BaseGroup;
        SuperGroup = g as ConcreteGroup<T>;
        if (singleton)
        {
            Elements = new HashSet<T> { ne };
            var kp = new Dictionary<T, int> { [ne] = 1 };
            ElementsOrders = new ReadOnlyDictionary<T, int>(kp);
            var lc = new Dictionary<T, ReadOnlyDictionary<T, int>> { [ne] = ElementsOrders };
            LongestCycles = new ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>>(lc);
            GroupType = GroupType.AbelianGroup;
        }
        else
        {
            Elements = SuperGroup?.GetElements().ToHashSet() ?? new HashSet<T>(g);
            LongestCycles = SuperGroup?.LongestCycles ?? Group.LongestCycles(g, Elements);
            ElementsOrders = SuperGroup?.ElementsOrders ?? Group.ElementsOrders(LongestCycles);
            GroupType = SuperGroup?.GroupType ?? (Group.IsCommutative(g, LongestCycles.Keys)
                ? GroupType.AbelianGroup
                : GroupType.NonAbelianGroup);
        }
    }

    public ConcreteGroup(IGroup<T> g) : this(g.Name, g)
    {
    }

    public ConcreteGroup(string name, IGroup<T> g, T[] generators)
    {
        Name = name;
        Hash = new HashCode().ToHashCode();
        var ne = g.Neutral();
        BaseGroup = ne.BaseGroup;
        SuperGroup = g as ConcreteGroup<T>;
        if (SuperGroup is not null && generators.Any(e => !SuperGroup.Contains(e)))
            throw new GroupException(GroupExceptionType.GroupDef);

        var tmpElements = Group.GenerateElements(g, generators);
        Elements = new HashSet<T>(tmpElements);
        LongestCycles = Group.LongestCycles(g, Elements);
        ElementsOrders = Group.ElementsOrders(LongestCycles);
        GroupType = Group.IsCommutative(g, LongestCycles.Keys)
            ? GroupType.AbelianGroup
            : GroupType.NonAbelianGroup;
    }

    public ConcreteGroup(IGroup<T> g, T[] generators) : this(g.Name, g, generators)
    {
    }

    protected HashSet<T> Elements { get; set; }
    public IGroup<T> BaseGroup { get; }
    public IConcreteGroup<T>? SuperGroup { get; }
    public GroupType GroupType { get; protected set; }
    public string Name { get; protected set; }
    public ReadOnlyDictionary<T, int> ElementsOrders { get; protected set; }
    public ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>> LongestCycles { get; protected set; }

    public IEnumerator<T> GetEnumerator()
    {
        return GetElements().GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetElements().GetEnumerator();
    }

    public bool Equals(IGroup<T>? other)
    {
        return other?.Hash == Hash;
    }

    public int Hash { get; }

    public IEnumerable<T> GetElements()
    {
        return Elements;
    }

    public virtual T Neutral()
    {
        return SuperGroup?.Neutral() ?? BaseGroup.Neutral();
    }

    public virtual T Invert(T e)
    {
        return SuperGroup?.Invert(e) ?? BaseGroup.Invert(e);
    }

    public virtual T Op(T e1, T e2)
    {
        return SuperGroup?.Op(e1, e2) ?? BaseGroup.Op(e1, e2);
    }

    public T this[params ValueType[] us]
    {
        get
        {
            var e = BaseGroup[us];
            if (!this.Contains(e))
                throw new GroupException(GroupExceptionType.GroupDef);

            return e;
        }
    }

    public bool IsIsomorphicTo<TU>(IConcreteGroup<TU> gu) where TU : struct, IElt<TU>
    {
        var set0 = ElementsOrders.Values.Ascending();
        var set1 = gu.ElementsOrders.Values.Ascending();
        return set0.SequenceEqual(set1);
    }

    public override int GetHashCode()
    {
        return Hash;
    }

    public override string ToString()
    {
        return Name;
    }
}