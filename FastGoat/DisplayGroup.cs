namespace FastGoat;

public enum SortBy
{
    Value,
    Order
}

public static class DisplayGroup
{
    public static void Head<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        Console.WriteLine($"|{g.Name}| = {g.Count()}");
        Console.WriteLine($"Type        {g.GroupType}");
        Console.WriteLine($"BaseGroup   {g.BaseGroup}");
        if (g.SuperGroup is not null)
        {
            var superGroup = g.SuperGroup;
            Console.WriteLine($"SuperGroup  |{superGroup}| = {superGroup.Count()}");
        }

        if (g is QuotientGroup<T> quo)
        {
            var normalGroup = quo.NormalSubGroup;
            Console.WriteLine($"NormalGroup |{normalGroup}| = {normalGroup.Count()}");
        }

        Console.WriteLine();
    }

    public static void Elements<T>(ConcreteGroup<T> g, SortBy sortBy = SortBy.Order) where T : struct, IElt<T>
    {
        if (g.ElementsOrders.Count > 300)
        {
            Console.WriteLine("   ***   Too Big   ***");
            return;
        }

        var ordered = sortBy == SortBy.Value
            ? g.ElementsOrders.Keys.Ascending()
            : g.ElementsOrders.Keys.OrderBy(a => g.ElementsOrders[a]).ThenAscending();

        var k = 0;
        var digits1 = $"{g.ElementsOrders.Count}".Length;
        var digits2 = $"{g.ElementsOrders.Values.Max()}".Length;
        var fmt = $"({{0,{digits1}}})[{{1,{digits2}}}] = {{2}}";

        Console.WriteLine("Elements");
        foreach (var elt in ordered)
            Console.WriteLine(fmt, ++k, g.ElementsOrders[elt], elt);

        Console.WriteLine();
    }

    public static void Classes<T>(QuotientGroup<T> g) where T : struct, IElt<T>
    {
        Console.WriteLine("Cosets");
        if (g.Representatives.Values.Count() > 40)
        {
            Console.WriteLine("   ***   Too Big   ***");
            return;
        }

        var ordered = g.ElementsOrders.Keys.OrderBy(a => g.ElementsOrders[a]).ThenAscending();

        var k = 0;
        var digits1 = $"{g.ElementsOrders.Count}".Length;
        var digits2 = $"{g.ElementsOrders.Values.Max()}".Length;
        var fmt = $"({{0,{digits1}}})[{{1,{digits2}}}]";

        foreach (var elt in ordered)
        {
            Console.WriteLine(fmt, ++k, g.ElementsOrders[elt]);
        }

        Console.WriteLine();
    }

    public static void Cosets<T>(QuotientGroup<T> g) where T : struct, IElt<T>
    {
        Console.WriteLine("Cosets");
        if (g.Representatives.Values.Count() > 40)
        {
            Console.WriteLine("   ***   Too Big   ***");
            return;
        }

        var ordered = g.ElementsOrders.Keys.OrderBy(a => g.ElementsOrders[a]).ThenAscending();

        var k = 0;
        var digits1 = $"{g.ElementsOrders.Count}".Length;
        var digits2 = $"{g.ElementsOrders.Values.Max()}".Length;
        var fmt = $"({{0,{digits1}}})[{{1,{digits2}}}]";

        foreach (var elt in ordered)
        {
            Console.WriteLine(fmt, ++k, g.ElementsOrders[elt]);
            var set = g.Cosets[elt];
            foreach (var elt1 in set.Ascending()) Console.WriteLine($"    {elt1}");
        }

        Console.WriteLine();
    }

    public static void Table<T>(ConcreteGroup<T> g, SortBy sortBy = SortBy.Order) where T : struct, IElt<T>
    {
        if (g.ElementsOrders.Count > 300)
        {
            Console.WriteLine("   ***   Too Big   ***");
            return;
        }

        var ordered = sortBy == SortBy.Value
            ? g.ElementsOrders.Keys.Ascending().ToArray()
            : g.ElementsOrders.Keys.OrderBy(a => g.ElementsOrders[a]).ThenAscending().ToArray();

        var digits1 = $"{g.ElementsOrders.Count}".Length;
        var fmt = $"{{0,{digits1}}}";
        var elt2Int = ordered.Select((e, i) => (e, i)).ToDictionary(p => p.e, p => p.i + 1);
        var rows = ordered.Select(e1 => ordered.Select(e2 => elt2Int[g.Op(e1, e2)]));
        foreach (var row in rows)
        {
            Console.WriteLine(row.Glue(" ", fmt));
        }

        Console.WriteLine();
    }

    public static void HeadElements<T>(ConcreteGroup<T> g, SortBy sortBy = SortBy.Order) where T : struct, IElt<T>
    {
        Head(g);
        Elements(g, sortBy);
    }

    public static void HeadTable<T>(ConcreteGroup<T> g, SortBy sortBy = SortBy.Order) where T : struct, IElt<T>
    {
        Head(g);
        Table(g, sortBy);
    }

    public static void HeadElementsTable<T>(ConcreteGroup<T> g, SortBy sortBy = SortBy.Order) where T : struct, IElt<T>
    {
        Head(g);
        Elements(g, sortBy);
        Table(g, sortBy);
    }

    public static void HeadClasses<T>(QuotientGroup<T> g) where T : struct, IElt<T>
    {
        Head(g);
        Classes(g);
    }

    public static void HeadCosets<T>(QuotientGroup<T> g) where T : struct, IElt<T>
    {
        Head(g);
        Cosets(g);
    }

    public static void HeadClassesTable<T>(QuotientGroup<T> g) where T : struct, IElt<T>
    {
        Head(g);
        Classes(g);
        Table(g, SortBy.Order);
    }

    public static void HeadCosetsTable<T>(QuotientGroup<T> g) where T : struct, IElt<T>
    {
        Head(g);
        Cosets(g);
        Table(g, SortBy.Order);
    }

    public static void HeadSdp<T1, T2>(SemiDirectProduct<T1, T2> p)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        Console.WriteLine($"|{p.Name}| = {p.Count()}");
        Console.WriteLine($"Type        {p.GroupType}");
        Console.WriteLine($"BaseGroup    {p.BaseGroup}");
        Console.WriteLine($"NormalGroup  |{p.N}| = {p.N.Count()}");
        Console.WriteLine($"ActionGroup  |{p.G}| = {p.G.Count()}");

        Console.WriteLine();
        Console.WriteLine("Actions");

        var ordered = p.G.ElementsOrders.Keys.OrderBy(a => p.G.ElementsOrders[a]).ThenAscending().ToArray();
        foreach (var g in ordered)
        {
            Console.WriteLine("g={0} y(g):{1}", g, p.ActionsStr[g]);
        }

        Console.WriteLine();
    }

    public static void HeadElementsSdp<T1, T2>(SemiDirectProduct<T1, T2> p, SortBy sortBy = SortBy.Order)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        HeadSdp(p);
        Elements(p, sortBy);
    }

    public static void HeadTableSdp<T1, T2>(SemiDirectProduct<T1, T2> p, SortBy sortBy = SortBy.Order)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        HeadSdp(p);
        Table(p, sortBy);
    }

    public static void HeadElementsTableSdp<T1, T2>(SemiDirectProduct<T1, T2> p, SortBy sortBy = SortBy.Order)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        HeadSdp(p);
        Elements(p, sortBy);
        Table(p, sortBy);
    }
}