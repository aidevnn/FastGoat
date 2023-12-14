namespace FastGoat.Structures.GenericGroup;

public class IsomorphEquality<T> : EqualityComparer<ConcreteGroup<T>> where T : struct, IElt<T>
{
    public override bool Equals(ConcreteGroup<T>? x, ConcreteGroup<T>? y)
    {
        return x is not null && y is not null && x.IsIsomorphicTo(y);
    }

    public override int GetHashCode(ConcreteGroup<T> obj) => (obj.Count(), obj.GroupType).GetHashCode();
}

public class IsomorphSubGroupsInfosEquality<T> : EqualityComparer<AllSubgroups<T>> where T : struct, IElt<T>
{
    public override bool Equals(AllSubgroups<T> x, AllSubgroups<T> y)
    {
        return x.Parent.IsIsomorphicTo(y.Parent);
    }

    public override int GetHashCode(AllSubgroups<T> obj) => obj.Infos.GetHashCode();
}