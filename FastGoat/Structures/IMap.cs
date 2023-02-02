using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures;

public interface IMap<T1, T2> : IElt<IMap<T1, T2>> where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
{
    ConcreteGroup<T1> Domain { get; }
    IEnumerable<T1> Kernel();
    IEnumerable<T2> Image();
    T2 this[T1 index] { get; }

    public static bool EqualiltyMap(IMap<T1, T2> map1, IMap<T1, T2> map2)
    {
        if (!map1.Domain.SetEquals(map2.Domain))
            return false;

        return map1.Domain.All(e => map1[e].Equals(map2[e]));
    }
    
    public static int CompareMap(IMap<T1, T2> map1, IMap<T1, T2> map2)
    {
        if (!map1.Domain.SetEquals(map2.Domain))
            return 1;

        foreach (var e in map1.Domain.Ascending())
        {
            var compE = map1[e].CompareTo(map2[e]);
            if (compE != 0)
                return compE;
        }

        return 0;
    }
}