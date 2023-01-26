using System.Collections;
using System.Collections.ObjectModel;
using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public class ConcreteGroup<T> : IGroup<T> where T : struct, IElt<T>
{
    public ConcreteGroup(string name, IGroup<T> g, bool singleton = false)
    {
        Name = name;
        Hash = Guid.NewGuid().GetHashCode();
        var ne = g.Neutral();
        SuperGroup = g as ConcreteGroup<T>;
        BaseGroup = SuperGroup?.BaseGroup ?? g;
        if (singleton)
        {
            Elements = new HashSet<T> { ne };
            PseudoGenerators = new(new List<T>());
            var kp = new Dictionary<T, int> { [ne] = 1 };
            ElementsOrders = new ReadOnlyDictionary<T, int>(kp);
            var lc = new Dictionary<T, ReadOnlyDictionary<T, int>> { [ne] = ElementsOrders };
            GroupType = GroupType.AbelianGroup;
        }
        else
        {
            var pseudo = SuperGroup?.GetGenerators().ToList() ?? g.GetGenerators().ToList();
            PseudoGenerators = new(pseudo);
            Elements = SuperGroup?.GetElements().ToHashSet() ?? new HashSet<T>(g);
            ElementsOrders = SuperGroup?.ElementsOrders ?? Group.ElementsOrders(g, Elements);
            GroupType = SuperGroup?.GroupType ?? (Group.IsCommutative(g, PseudoGenerators)
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
        Hash = Guid.NewGuid().GetHashCode();
        var ne = g.Neutral();
        SuperGroup = g as ConcreteGroup<T>;
        BaseGroup = SuperGroup?.BaseGroup ?? g;
        if (SuperGroup is not null && generators.Any(e => !SuperGroup.Contains(e)))
            throw new GroupException(GroupExceptionType.GroupDef);

        var (tmpElements, uniqueGenerators) = Group.UniqueGenerators(this, generators);
        Elements = new HashSet<T>(tmpElements);
        ElementsOrders = Group.ElementsOrders(g, Elements);
        GroupType = Group.IsCommutative(g, generators)
            ? GroupType.AbelianGroup
            : GroupType.NonAbelianGroup;

        PseudoGenerators = new(uniqueGenerators);
    }

    public ConcreteGroup(IGroup<T> g, T[] generators) : this(g.Name, g, generators)
    {
    }

    public string ShortName => $"|{Name}| = {this.Count()}";
    public ReadOnlyCollection<T> PseudoGenerators { get; protected set; }
    public IEnumerable<int> ElementsOrdersList() => ElementsOrders.Values.Ascending();
    protected HashSet<T> Elements { get; set; }
    public bool SetEquals(IEnumerable<T> ts) => Elements.SetEquals(ts);
    public bool SubSetOf(IEnumerable<T> ts) => Elements.IsSubsetOf(ts);
    public IGroup<T> BaseGroup { get; }
    public ConcreteGroup<T>? SuperGroup { get; }
    public GroupType GroupType { get; protected set; }
    public string Name { get; protected set; }
    public virtual string[] Details => Array.Empty<string>();
    public ReadOnlyDictionary<T, int> ElementsOrders { get; protected set; }

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

    public int Hash { get; protected set; }

    public IEnumerable<T> GetGenerators()
    {
        if (Elements.Count == 1)
            yield return Neutral();
        else if (ElementsOrders.Values.Max() == Elements.Count())
            yield return ElementsOrders.MaxBy(p => p.Value).Key;
        else if (PseudoGenerators.Count() != 0)
        {
            foreach (var e in PseudoGenerators)
                yield return e;
        }
        else
        {
            throw new GroupException(GroupExceptionType.GroupDef);
        }
    }

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

    public bool IsIsomorphicTo<Tu>(ConcreteGroup<Tu> gu) where Tu : struct, IElt<Tu>
    {
        if (!ElementsOrdersList().Ascending().SequenceEqual(gu.ElementsOrdersList().Ascending()))
            return false;

        var isos = Group.AllMorphisms(this, gu, Group.MorphismType.Isomorphism);
        return isos.Any(h => h.Count != 0);
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