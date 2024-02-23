namespace FastGoat.Structures.Subgroups;

public record SubGroupsInfos(int AllSubGr, int AllConjsCl, int AllNorms) : IComparable<SubGroupsInfos>
{
    public (int, int, int) ToTuples() => (AllSubGr, AllConjsCl, AllNorms);

    public int CompareTo(SubGroupsInfos? other)
    {
        if (other is null)
            return 1;

        return ToTuples().CompareTo(other.ToTuples());
    }
}