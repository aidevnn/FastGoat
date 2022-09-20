namespace FastGoat;

public partial class ConcreteGroup<T>
{
    protected enum SubGroupEnum { OutOfBaseGroup, OutOfSubGroup, IsSubGroup }
    protected struct Order : IEquatable<Order>, IComparable<Order>
    {
        public T e { get; }
        public T g { get; }
        public int p { get; }
        public Order(T e0)
        {
            e = g = e0;
            p = 1;
        }
        public Order(T e0, T g0, int p0)
        {
            e = e0;
            g = g0;
            p = p0;
        }

        public int CompareTo(Order other) => p.CompareTo(other.p);
        public bool Equals(Order other) => e.Equals(other.e);
        public override int GetHashCode() => e.Hash;
        public override string ToString() => $"{e} = g*{p,-2}";

        public static implicit operator T(Order o) => o.e;
    }
    protected IEnumerable<T> Generate(IEnumerable<T> set)
    {
        return Generate(new[] { this.Neutral() }, set);
    }
    protected IEnumerable<T> Generate(IEnumerable<T> pset, IEnumerable<T> set)
    {
        // int count = 0;
        var n = this.Neutral();
        if (set.Count() == 0)
            return new[] { this.Neutral() };

        var set0 = set.ToHashSet();
        Queue<T> newElements = new(pset);
        HashSet<T> all = new(pset);

        while (newElements.Count != 0)
        {
            var e0 = newElements.Dequeue();
            foreach (var e1 in set0)
            {
                // ++count;
                var e2 = this.Op(e0, e1);
                if (all.Add(e2))
                    newElements.Enqueue(e2);
            }
        }

        // Console.WriteLine($"Counter : {count}; Total ; {all.Count}");
        return all.Ascending();
    }
    protected GroupType ComputeGroupType(IEnumerable<Order> gens)
    {
        foreach (var e0 in gens)
        {
            foreach (var e1 in gens)
            {
                var e01 = this.Op(e0.e, e1.e);
                var e10 = this.Op(e1.e, e0.e);
                if (!e01.Equals(e10))
                    return GroupType.NotAbelianGroup;
            }
        }
        return GroupType.AbelianGroup;
    }
    protected HashSet<Order> Monogenic(T e)
    {
        var n = this.Neutral();
        var e0 = e;
        var p = 1;
        HashSet<Order> set = new() { new(e) };
        while (!n.Equals(e0))
        {
            e0 = this.Op(e0, e);
            ++p;
            set.Add(new(e0, e, p));
        }

        return set;
    }
    protected Dictionary<T, int> ComputeOrders(IEnumerable<HashSet<Order>> gens)
    {
        HashSet<(T, int)> ordersTuple = new();
        foreach (var p in gens)
        {
            int m = p.Count;
            ordersTuple.UnionWith(p.Select(t => (t.e, m / IntExt.GCD(m, t.p))));
        }
        return ordersTuple.ToDictionary(a => a.Item1, a => a.Item2);
    }
    Dictionary<Order, HashSet<Order>> monogenicSubGroup;
    public Dictionary<T, (T, int)[]> MonogenicSubGroup() => monogenicSubGroup.ToDictionary(a => a.Key.e, b => b.Value.Select(e => (e.e, e.p)).Ascending().ToArray());
    protected Dictionary<Order, HashSet<Order>> ComputeGenerators(IEnumerable<T> elements)
    {
        var ne = this.Neutral();
        HashSet<Order> set = new(elements.Select(e => new Order(e)));
        Dictionary<Order, HashSet<Order>> gens = new();

        while (set.Count != 0)
        {
            var e = set.First();
            var g = Monogenic(e.e);
            set.ExceptWith(g);
            if (gens.Count == 0)
            {
                gens[e] = g;
                continue;
            }

            var gens0 = new Dictionary<Order, HashSet<Order>>(gens);
            gens.Clear();
            bool found = false;
            foreach (var p in gens0)
            {
                var e0 = p.Key;
                var g0 = p.Value;

                if (g.Count > g0.Count)
                {
                    if (g.IsSupersetOf(g0))
                    {
                        if (!gens.ContainsKey(e))
                        {
                            gens[e] = g;
                            found = true;
                        }
                    }
                    else
                        gens[e0] = g0;
                }
                else
                {
                    gens[e0] = g0;
                    if (!found && g0.IsSubsetOf(g))
                        found = true;
                }
            }

            if (!found)
                gens[e] = g;
        }

        return gens;
    }
    protected (SubGroupEnum, Dictionary<Order, HashSet<Order>>) VerifySubGroup(IEnumerable<T> ts)
    {
        if (ts.Any(t => !elements.Contains(t)))
            return (SubGroupEnum.OutOfBaseGroup, new());

        var gens = ComputeGenerators(ts);
        var ts0 = gens.SelectMany(p => p.Value).Select(p => p.e).ToHashSet();
        if (!ts0.SetEquals(ts))
            return (SubGroupEnum.OutOfSubGroup, new());

        var head = gens.Keys.Select(p => p.e).Ascending().ToHashSet();
        foreach (var e0 in head)
        {
            var ei = Invert(e0);
            if (!ts0.Contains(ei))
                return (SubGroupEnum.OutOfSubGroup, new());

            foreach (var e1 in head)
            {
                if (!ts0.Contains(this.Op(ei, e1)))
                    return (SubGroupEnum.OutOfSubGroup, new());
            }
        }
        return (SubGroupEnum.IsSubGroup, gens);
    }
    protected (GroupType, Dictionary<T, int>, List<T>) ComputeDetails(IEnumerable<T> elements)
    {
        var gens = ComputeGenerators(elements);
        return ComputeDetails(gens);
    }
    protected (GroupType, Dictionary<T, int>, List<T>) ComputeDetails(Dictionary<Order, HashSet<Order>> gens)
    {
        monogenicSubGroup = new(gens);
        var gType = ComputeGroupType(gens.Keys);
        var orders = ComputeOrders(gens.Values);
        var monogenics = gens.Select(e => e.Key.e).Ascending().ToList();

        return (gType, orders, monogenics);
    }
    public Dictionary<T, HashSet<T>> LeftCosets(IEnumerable<T> H)
    {
        return elements.ToDictionary(x => x, x => H.Select(h => this.Op(x, h)).ToHashSet());
    }
    public Dictionary<T, HashSet<T>> RightCosets(IEnumerable<T> H)
    {
        return elements.ToDictionary(x => x, x => H.Select(h => this.Op(h, x)).ToHashSet());
    }
    public IEnumerable<IEnumerable<T>> Cosets1(ConcreteGroup<T> H)
    {
        if (!H.IsProperSubGroupOf(this))
            throw new SubGroupException("Not subgroup");

        var leftCosets = this.LeftCosets(H);
        var rightCosets = this.RightCosets(H);
        if (this.Any(x => !leftCosets[x].SetEquals(rightCosets[x])))
            throw new SubGroupException("Not normal subgroup");

        var cosets = leftCosets.Values.Distinct(new EqualityHashSet<T>()).Select(set => set.Ascending());
        return cosets;
    }
    public IEnumerable<IEnumerable<T>> Cosets(ConcreteGroup<T> H)
    {
        if (!H.IsSubGroupOf(this))
            throw new SubGroupException("Not subgroup");

        List<HashSet<T>> cosets = new();
        foreach (var x in elements)
        {
            var xi = this.Invert(x);
            var xH = H.Select(h => this.Op(x, h)).ToHashSet();
            var xHxi = xH.Select(xh => this.Op(xh, xi));
            if (!H.GroupEqual(xHxi))
                throw new SubGroupException("Not normal subgroup");

            cosets.Add(xH);
        }

        var fcosets = cosets.Distinct(new EqualityHashSet<T>()).Select(set => set.Ascending());
        return fcosets;
    }
}
