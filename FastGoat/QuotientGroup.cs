namespace FastGoat;

public class QuotientGroup<T> : WorkGroup<T> where T : struct, IElt<T>
{
    public QuotientGroup(WorkGroup<T> grG, WorkGroup<T> grH) : base(grG)
    {
        if (!grH.IsSubGroupOf(grG))
        {
            grG.DisplayHead("G");
            grH.DisplayHead("H");
            grH.DisplayElements(SortElements.ByValue);
            throw new Exception("H is not a subgroup of G");
        }

        G = new(grG);
        InternalTable = Table(grH);

        var tmpElements = G.SortByOrder(InternalTable.Values.ToHashSet());
        (groupType, elementOrder) = ComputeDetails(tmpElements);
        elements = new(tmpElements);
    }
    WorkGroup<T> G { get; }
    Dictionary<T, T> InternalTable { get; }
    Dictionary<T, T> Table(WorkGroup<T> H)
    {
        Dictionary<T, T> table = new();

        var cosets = G.Cosets(H).Select(h => h.Ascending());
        foreach (var set in cosets)
        {
            var r = set.First();
            foreach (var e in set)
                table[e] = r;
        }

        return table;
    }

    public override T Neutral() => G.Neutral();
    public override T Invert(T a) => InternalTable[G.Invert(a)];
    public override T Op(T a, T b) => InternalTable[G.Op(a, b)];

    public void DisplayCosets(SortElements sort = SortElements.ByOrder)
    {
        var cosets = InternalTable.GroupBy(a => a.Value).ToDictionary(a => a.Key, a => a.Select(b => b.Key).Ascending());
        int k = 0;
        var arr = sort == SortElements.ByOrder ? SortByOrder(elements).ToArray() : elements.Ascending().ToArray();
        foreach (var a in arr)
        {
            DisplayElement(a, $"{++k}");
            foreach (var e in cosets[a])
                Console.WriteLine($"      {e}");
        }
    }
}
