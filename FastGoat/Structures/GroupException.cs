namespace FastGoat.Structures;

public enum GroupExceptionType
{
    GroupDef,
    BaseGroup,
    NotSubGroup,
    NotNormal,
    OnlyCyclicGroups,
    SemiDirectProductDontExist
}

public class GroupException : Exception
{
    public GroupException(GroupExceptionType type) : base(ErrMsg(type))
    {
    }

    private static string ErrMsg(GroupExceptionType type)
    {
        return type switch
        {
            GroupExceptionType.GroupDef => "Group or Element cannot be created",
            GroupExceptionType.BaseGroup => "Groups or Elements do not belong to the base group",
            GroupExceptionType.NotSubGroup => "Not a valid subgroup",
            GroupExceptionType.NotNormal => "Not a normal subgroup",
            GroupExceptionType.OnlyCyclicGroups =>
                "Only cyclic groups are allowed at this time for semi-direct product",
            GroupExceptionType.SemiDirectProductDontExist => "Semi-direct product does not exists",
            _ => "Unexpected error"
        };
    }
}