using System.Collections;
using System.Collections.ObjectModel;

namespace FastGoat;

public class ConcreteGroup<T> : IConcreteGroup<T> where T : struct, IElt<T>
{
    public ConcreteGroup(string name, IGroup<T> g, bool singleton = false)
    {
        Name = name;
        Hash = Guid.NewGuid().GetHashCode();
        var ne = g.Neutral();
        BaseGroup = ne.BaseGroup;
        SuperGroup = g as ConcreteGroup<T>;
        if (singleton)
        {
            Elements = new HashSet<T> { ne };
            PseudoGenerators = new List<T>();
            var kp = new Dictionary<T, int> { [ne] = 1 };
            ElementsOrders = new ReadOnlyDictionary<T, int>(kp);
            var lc = new Dictionary<T, ReadOnlyDictionary<T, int>> { [ne] = ElementsOrders };
            LongestCycles = new ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>>(lc);
            GroupType = GroupType.AbelianGroup;
        }
        else
        {
            PseudoGenerators = SuperGroup?.GetGenerators().ToList() ?? g.GetGenerators().ToList();
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
        Hash = Guid.NewGuid().GetHashCode();
        var ne = g.Neutral();
        BaseGroup = ne.BaseGroup;
        SuperGroup = g as ConcreteGroup<T>;
        if (SuperGroup is not null && generators.Any(e => !SuperGroup.Contains(e)))
            throw new GroupException(GroupExceptionType.GroupDef);

        var (tmpElements, uniqueGenerators) = InternalGenerators(generators);
        PseudoGenerators = uniqueGenerators;
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

    protected List<T> PseudoGenerators { get; set; }

    protected (HashSet<T> elements, List<T> uniqueGenerators) InternalGenerators(T[] generators)
    {
        HashSet<T> tmpElements = new() { this.Neutral() };
        List<T> uniqueGenerators = new();
        foreach (var elt in generators)
        {
            if (tmpElements.Contains(elt))
                continue;

            uniqueGenerators.Add(elt);
            tmpElements = Group.GenerateElements(this, tmpElements, uniqueGenerators);
        }

        return (tmpElements, uniqueGenerators);
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

    public int Hash { get; protected set; }

    public IEnumerable<T> GetGenerators()
    {
        if (LongestCycles.Count == 1)
            yield return LongestCycles.First().Key;
        else if (PseudoGenerators.Count() != 0)
        {
            foreach (var e in PseudoGenerators)
                yield return e;
        }
        else
        {
            var bgGens = BaseGroup.GetGenerators().ToArray();
            var prod = bgGens.Aggregate(1, (acc, a) => ElementsOrders[a] * acc);
            if (prod == Elements.Count)
            {
                foreach (var e in bgGens)
                    yield return e;
            }
            else
            {
                throw new GroupException(GroupExceptionType.GroupDef);
            }
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
        var homs = Group.AllHomomorphisms(this, gu);
        var nbIsomorphisms = homs.Count(h => h.Values.Distinct().Count() == h.Count);
        return nbIsomorphisms > 0;
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