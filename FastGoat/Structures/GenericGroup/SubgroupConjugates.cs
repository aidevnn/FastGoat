using System.Drawing;

namespace FastGoat.Structures.GenericGroup;

public readonly struct SubgroupConjugates<T> : IElt<SubgroupConjugates<T>> where T : struct, IElt<T>
{
    public ConcreteGroup<T> Parent { get; }
    public List<ConcreteGroup<T>> Conjugates { get; }
    public ConcreteGroup<T> Representative { get; }
    public int Index { get; }
    public int Order { get; }

    public SubgroupConjugates(ConcreteGroup<T> parent, ConcreteGroup<T> subGroup)
    {
        Parent = parent;
        Conjugates = Group.SubGroupsConjugates(parent, subGroup);
        Representative = subGroup;
        Order = Representative.Count();
        Index = parent.Count() / Order;
        Hash = (Parent.Hash, Order, Index).GetHashCode();
        if (Order == 1)
            Representative.Name = "()";
    }

    public SubgroupConjugates(ConcreteGroup<T> parent, List<ConcreteGroup<T>> conjugates)
    {
        Parent = parent;
        Conjugates = conjugates.Count == 0 ? new() { Group.Generate("()", parent, parent.Neutral()) } : conjugates;
        Representative = Conjugates[0];
        Order = Representative.Count();
        Index = parent.Count() / Order;
        Hash = (Parent.Hash, Order, Index).GetHashCode();
    }

    public SubgroupConjugates<T> Restriction(ConcreteGroup<T> g) => new(g, Conjugates.Where(e => e.SubSetOf(g)).ToList());
    public bool Contains(ConcreteGroup<T> g) => Conjugates.Any(e => e.SetEquals(g));

    public GroupType GroupType => Representative.GroupType;
    public int Size => Conjugates.Count;
    public (int, int, GroupType) OST => (Order, Size, GroupType);
    public bool IsNormal => Size == 1;
    public bool IsProperNormal => IsNormal && Order != 1 && Index != 1;
    public bool Equals(SubgroupConjugates<T> other) => Hash == other.Hash && Conjugates.Any(e => e.SetEquals(other.Representative));

    public int CompareTo(SubgroupConjugates<T> other) => OST.CompareTo(other.OST);

    public int Hash { get; }

    public override int GetHashCode() => Hash;
    public override string ToString() => Representative.Name;
}