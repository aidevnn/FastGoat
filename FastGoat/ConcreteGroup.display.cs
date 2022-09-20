namespace FastGoat;

public partial class ConcreteGroup<T>
{
    public virtual void DisplayHead(string name = "G")
    {
        Console.WriteLine("|{0}| = {1} in {2}", name, this.Count(), BaseGroup);
        Console.WriteLine($"is {groupType}");
    }

    public void DisplayElement(T e, string symb = "@", string fmt = "({0})[{1}] = {2}") => Console.WriteLine(fmt, symb, GetOrderOf(e), e);

    public void DisplayElements(SortElements sort = SortElements.ByOrder)
    {
        Console.WriteLine("Elements");
        if (elements.Count > 300)
        {
            Console.WriteLine("*** TOO BIG ***");
            return;
        }

        var arr = sort == SortElements.ByOrder ? SortByOrder(elements).ToArray() : elements.Ascending().ToArray();
        var digits = arr.Length.ToString().Length;
        var fmt = $"({{0,{digits}}})[{{1,{digits}}}] = {{2}}";
        for (int k = 0; k < arr.Length; ++k)
        {
            var e = arr[k];
            DisplayElement(e, $"{k + 1}", fmt);
            // var o = GetOrderOf(e);
            // Console.WriteLine(fmt, k + 1, o, e);
        }
    }

    public void DisplayTable(SortElements sort = SortElements.ByOrder)
    {
        Console.WriteLine("Table");
        if (elements.Count > 40)
        {
            Console.WriteLine("*** TOO BIG ***");
            return;
        }

        var arr = sort == SortElements.ByOrder ? SortByOrder(elements).ToArray() : elements.Ascending().ToArray();
        var dico = arr.Select((e, i) => (e, i + 1)).ToDictionary(a => a.e, a => a.Item2);
        var table = arr.Select(a0 => arr.Select(a1 => dico[this.Op(a0, a1)]).ToArray()).ToArray();
        var digits = arr.Length.ToString().Length;
        var fmt = $"{{0,{digits}}}";
        foreach (var r in table)
            Console.WriteLine(r.Glue(" ", fmt));
    }
    public (int[][], int[][]) CayleyTable(SortElements sort = SortElements.ByOrder)
    {
        var arr = sort == SortElements.ByOrder ? SortByOrder(elements).ToArray() : elements.Ascending().ToArray();
        var dico = arr.Select((e, i) => (e, i + 1)).ToDictionary(a => a.e, a => a.Item2);
        var rows = arr.Select(a0 => arr.Select(a1 => dico[this.Op(a0, a1)]).ToArray()).ToArray();
        var cols = arr.Select(a1 => arr.Select(a0 => dico[this.Op(a0, a1)]).ToArray()).ToArray();
        return (rows, cols);
    }
    public virtual void DisplayDetails(SortElements sort) => DisplayDetails("G", sort);
    public virtual void DisplayDetails(string name = "G", SortElements sort = SortElements.ByOrder)
    {
        DisplayHead(name);
        Console.WriteLine();
        DisplayElements(sort);
        Console.WriteLine();
        DisplayTable(sort);
        Console.WriteLine();
    }
}
