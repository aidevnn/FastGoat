namespace FastGoat.Structures.GenericGroup;

public class IsomorphEquality<T> : EqualityComparer<ConcreteGroup<T>> where T : struct, IElt<T>
{
    public override bool Equals(ConcreteGroup<T>? x, ConcreteGroup<T>? y)
    {
        return x is not null && y is not null && x.IsIsomorphicTo(y);
    }

    public override int GetHashCode(ConcreteGroup<T> obj) => (obj.Count(), obj.GroupType).GetHashCode();
}