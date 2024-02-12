using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;

namespace FastGoat.Structures;

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

        Console.WriteLine();
    }

    public static void Head<T>(ConcreteGroup<Coset<T>> g) where T : struct, IElt<T>
    {
        Console.WriteLine($"|{g.Name}| = {g.Count()}");
        Console.WriteLine($"Type        {g.GroupType}");
        Console.WriteLine($"BaseGroup   {g.BaseGroup}");

        var g0 = (Quotient<T>)g.BaseGroup;
        Console.WriteLine($"Group           |{g0.G}| = {g0.G.Count()}");
        Console.WriteLine($"NormalSubGroup  |{g0.H}| = {g0.H.Count()}");

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
        var fmt = $"({{0,{digits1}}})[{{1,{digits2}}}] = ";

        Console.WriteLine("Elements");
        foreach (var elt in ordered)
        {
            Console.Write(fmt, ++k, g.ElementsOrders[elt]);
            Console.WriteLine(elt);
        }

        Console.WriteLine();
    }

    public static void Cosets<T>(ConcreteGroup<Coset<T>> g, bool details = false) where T : struct, IElt<T>
    {
        Console.WriteLine("Cosets");
        var g0 = (Quotient<T>)g.BaseGroup;
        if (g0.G.Count() > 40)
        {
            Console.WriteLine("   ***   Too Big   ***");
            return;
        }

        var ordered = g.ElementsOrders.Keys.OrderBy(a => g.ElementsOrders[a]).ThenAscending();

        var k = 0;
        var digits1 = $"{g.ElementsOrders.Count}".Length;
        var digits2 = $"{g.ElementsOrders.Values.Max()}".Length;
        var fmt = $"({{0,{digits1}}})[{{1,{digits2}}}] = {{2}}";

        foreach (var elt in ordered)
        {
            Console.WriteLine(fmt, ++k, g.ElementsOrders[elt], elt);
            if (details)
            {
                foreach (var elt1 in elt)
                    Console.WriteLine($"    {elt1}");
            }
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

    public static void Orders<T>(ConcreteGroup<T> g) where T : struct, IElt<T>
    {
        Console.WriteLine("Elements Orders : {0}",
            g.ElementsOrdersList().GroupBy(a => a).ToDictionary(a => a.Key, a => a.Count()).AscendingByKey()
                .GlueMap(fmt: "[{0}]:{1}"));
    }

    public static void HeadOrders<T>(ConcreteGroup<T> g, bool newline = true) where T : struct, IElt<T>
    {
        Head(g);
        Orders(g);
        if (newline)
            Console.WriteLine();
    }

    public static void HeadOrders<T>(ConcreteGroup<Coset<T>> g, bool newline = true) where T : struct, IElt<T>
    {
        Head(g);
        Orders(g);
        if (newline)
            Console.WriteLine();
    }

    public static void HeadSdpOrders<T1, T2>(SemiDirectProduct<T1, T2> g, bool newline = true)
        where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        HeadSdp(g);
        Orders(g);
        if (newline)
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

    public static void HeadCosets<T>(ConcreteGroup<Coset<T>> g, bool details = false) where T : struct, IElt<T>
    {
        Head(g);
        Cosets(g, details);
    }

    public static void HeadCosetsTable<T>(ConcreteGroup<Coset<T>> g) where T : struct, IElt<T>
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
        var faithfull = p.IsFaithFull() ? "FaithFull" : "Not FaithFull";
        Console.WriteLine($"Action {faithfull}");

        var ordered = p.G.ElementsOrders.Keys.Ascending().ThenBy(a => p.G.ElementsOrders[a]).ToArray();
        foreach (var g in ordered)
        {
            Console.WriteLine("g={0} y(g) = ({1})", g, p.Theta[g]);
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

    public static void AreIsomorphics<T1, T2>(ConcreteGroup<T1> g1, ConcreteGroup<T2> g2) where T1 : struct, IElt<T1>
        where T2 : struct, IElt<T2>
    {
        Console.WriteLine("{0} IsIsomorphicTo {1} : {2}", g1, g2, g1.IsIsomorphicTo(g2));
    }

    public static void ConjugacyClasses<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        Head(gr);
        new ConjugacyClasses<T>(gr).Display();
    }

    public static void HeadNames<T>(ConcreteGroup<T> g, bool setName = true) where T : struct, IElt<T>
    {
        var subGroups = new AllSubgroups<WElt>(g.ToGroupWrapper());
        var names = NamesTree.BuildName(subGroups);
        if (setName)
            g.Name = names[0].Name;

        HeadNames(g, subGroups.Infos, names);
    }

    public static void HeadElementsNames<T>(ConcreteGroup<T> g, bool setName = true) where T : struct, IElt<T>
    {
        var subGroups = new AllSubgroups<WElt>(g.ToGroupWrapper());
        var names = NamesTree.BuildName(subGroups);
        if (setName)
            g.Name = names[0].Name;

        HeadElementsNames(g, subGroups.Infos, names);
    }

    public static void HeadOrdersNames<T>(ConcreteGroup<T> g, bool setName = true) where T : struct, IElt<T>
    {
        var subGroups = new AllSubgroups<WElt>(g.ToGroupWrapper());
        var names = NamesTree.BuildName(subGroups);
        if (setName)
            g.Name = names[0].Name;

        HeadOrdersNames(g, subGroups.Infos, names);
    }

    public static void HeadNames<T>(ConcreteGroup<T> g, SubGroupsInfos infos, ANameElt[] names) where T : struct, IElt<T>
    {
        Head(g);
        Console.WriteLine(infos);
        names.Println("Group names");
        Console.WriteLine();
    }

    public static void HeadElementsNames<T>(ConcreteGroup<T> g, SubGroupsInfos infos, ANameElt[] names)
        where T : struct, IElt<T>
    {
        HeadElements(g);
        Console.WriteLine(infos);
        names.Println("Group names");
        Console.WriteLine();
    }

    public static void HeadOrdersNames<T>(ConcreteGroup<T> g, SubGroupsInfos infos, ANameElt[] names)
        where T : struct, IElt<T>
    {
        HeadOrders(g, newline: false);
        Console.WriteLine(infos);
        names.Println("Group names");
        Console.WriteLine();
    }
}