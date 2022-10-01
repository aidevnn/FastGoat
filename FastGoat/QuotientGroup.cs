using System.Collections.ObjectModel;

namespace FastGoat;

public class QuotientGroup<T> : ConcreteGroup<T> where T : struct, IElt<T>
{
    public QuotientGroup(ConcreteGroup<T> superGroup, ConcreteGroup<T> subGroup, string name) : base(name, superGroup,
        true)
    {
        NormalSubGroup = subGroup;
        var cosets = Group.Cosets(superGroup, subGroup);
        Cosets = new ReadOnlyDictionary<T, ReadOnlyCollection<T>>(cosets);
        var representatives =
            Cosets.SelectMany(p => p.Value.Select(v => (p.Key, v))).ToDictionary(p => p.v, p => p.Key);
        Representatives = new ReadOnlyDictionary<T, T>(representatives);
        Elements = Representatives.Values.ToHashSet();
        LongestCycles = Group.LongestCycles(this, Elements);
        ElementsOrders = Group.ElementsOrders(LongestCycles);
    }

    public ReadOnlyDictionary<T, ReadOnlyCollection<T>> Cosets { get; }
    public ReadOnlyDictionary<T, T> Representatives { get; }
    public ConcreteGroup<T> NormalSubGroup { get; }

    public override T Invert(T e)
    {
        var ei = base.Invert(e);
        return Representatives[ei];
    }

    public override T Op(T e1, T e2)
    {
        var e3 = base.Op(e1, e2);
        return Representatives[e3];
    }
}