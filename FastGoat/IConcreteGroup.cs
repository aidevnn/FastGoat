using System.Collections.ObjectModel;

namespace FastGoat;

public enum GroupType
{
    AbelianGroup,
    NonAbelianGroup
}

public interface IConcreteGroup<T> : IGroup<T> where T : struct, IElt<T>
{
    IGroup<T> BaseGroup { get; }
    IConcreteGroup<T>? SuperGroup { get; }
    ReadOnlyDictionary<T, int> ElementsOrders { get; }
    ReadOnlyDictionary<T, ReadOnlyDictionary<T, int>> LongestCycles { get; }
    GroupType GroupType { get; }
}