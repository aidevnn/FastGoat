namespace FastGoat;

public static class MapExt
{
    public static int CompareMapTo<T1, T2>(this IMap<T1, T2> map1, IMap<T1, T2> map2)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
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