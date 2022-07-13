using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;


public abstract class SubGroup<U> : SubSet<U>, ISubGroup<U> where U : struct, IElt<U>
{
    protected SubGroup(IGroup<U> group) : base(group)
    {
        UpperGroup = group;
        Infos = new DisplayGroup<U>(this);
        ElementOrder = new Dictionary<U, int>();
    }

    public IGroup<U> UpperGroup { get; }

    public abstract U Neutral { get; }
    public abstract U Invert(U a);
    public abstract U Op(U a, U b);

    Dictionary<U, int> ElementOrder { get; set; }
    public SortBy SortBy { get; set; } = SortBy.Order;

    public int GetOrder(U e) => ElementOrder[e];
    public bool EltOrder { get; private set; }
    public override void AddElement(U e)
    {
        EltOrder = false;
        base.AddElement(e);
    }

    public bool IsGroup()
    {
        foreach (var e0 in Elts)
        {
            var inv = Invert(e0);
            if (!Elts.Contains(inv))
                return false;
            foreach (var e1 in Elts)
            {
                if (!Elts.Contains(Op(e0, e1)))
                    return false;
            }
        }

        return true;
    }

    public bool IsCommutative()
    {
        foreach (var e0 in Elts)
        {
            foreach (var e1 in Elts)
            {
                var e2 = Op(e0, e1);
                if (!Elts.Contains(e2))
                    return false;

                var e3 = Op(e1, e0);
                if (!Elts.Contains(e3))
                    return false;

                if (!e2.Equals(e3))
                    return false;
            }
        }

        return true;
    }

    public bool DisplayOrders => ElementOrder.Count == 0 || ElementOrder.Values.All(a => a == 1);
    public void ComputeOrders()
    {
        ElementOrder = Elts.ToDictionary(a => a, b => 1);
        EltOrder = false;
        if (SortBy == SortBy.Value)
            return;

        List<(U e, int o)> orders = new List<(U, int)>();
        foreach (var e in Elts)
        {
            int ord = 0;
            var acc = Neutral;
            while (ord == 0 || !acc.Equals(Neutral))
            {
                ++ord;
                acc = Op(e, acc);
                if (!Elts.Contains(acc))
                    return;
            }

            orders.Add((e, ord));
        }

        EltOrder = true;
        ElementOrder = orders.ToDictionary(a => a.e, b => b.o);
    }

    public override int CompareElt(U a, U b)
    {
        if (SortBy == SortBy.Order && ElementOrder.ContainsKey(a) && ElementOrder.ContainsKey(b))
        {
            var ordA = ElementOrder[a];
            var ordB = ElementOrder[b];
            if (ordA != ordB)
                return ordA.CompareTo(ordB);
        }

        return base.CompareElt(a, b);
    }

    public void DisplayGroupTable(char symb = '*') => DisplayTable(symb, Op);

    public void Details(string? name = null, string? infos = null)
    {
        DisplayElements(name, infos);
        DisplayGroupTable();
    }
}