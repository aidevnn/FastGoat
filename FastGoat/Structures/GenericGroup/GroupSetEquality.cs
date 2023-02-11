namespace FastGoat.Structures.GenericGroup;

public class GroupSetEquality<T> : EqualityComparer<ConcreteGroup<T>> where T : struct, IElt<T>
{
    public override bool Equals(ConcreteGroup<T>? x, ConcreteGroup<T>? y)
    {
        return x is not null && y is not null && x.SetEquals(y);
    }

    public override int GetHashCode(ConcreteGroup<T> obj) => (obj.Count(), obj.GroupType).GetHashCode();
}