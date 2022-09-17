namespace FastGoat;

public partial class WorkGroup<T> : ConcreteGroup<T> where T : struct, IElt<T>
{
    public WorkGroup(IGroup<T> g) : base(g)
    {
        ControlGroup = new ConcreteGroup<T>(BaseGroup);
    }
    public WorkGroup(T e) : base(e.Group)
    {
        ControlGroup = new ConcreteGroup<T>(BaseGroup);
        var tmpElements = Generate(new[] { e });
        (groupType, elementOrder, monogenics) = ComputeDetails(tmpElements);
        elements = new(tmpElements);
    }
    public WorkGroup(IEnumerable<T> ts) : base(ts.First().Group)
    {
        ControlGroup = new ConcreteGroup<T>(BaseGroup);
        if (ts.Any(t => !t.Group.Equals(BaseGroup)))
            throw new BaseGroupException();

        var tmpElements = Generate(ts);
        (groupType, elementOrder, monogenics) = ComputeDetails(tmpElements);
        elements = new(tmpElements);
    }
    public WorkGroup(WorkGroup<T> group) : base(group.BaseGroup)
    {
        ControlGroup = group;
        (groupType, elementOrder) = (group.groupType, new(group.elementOrder));
        elements = group.elements.ToHashSet();
    }
    private WorkGroup(WorkGroup<T> group, IEnumerable<T> ts) : base(group.BaseGroup)
    {
        ControlGroup = group;
        if (ts.Any(e => !group.Contains(e)))
            throw new SubGroupException("Element doesnt belong to the super group");

        var tmpElements = Generate(ts);
        (groupType, elementOrder, monogenics) = ComputeDetails(tmpElements);
        elements = new(tmpElements);
    }
    private WorkGroup(WorkGroup<T> group, IEnumerable<T> ts, int diff) : base(group.BaseGroup)
    {
        ControlGroup = group;
        if (ts.Any(e => !group.BaseGroup.Equals(e.Group)))
            throw new BaseGroupException();

        var tmpElements = Generate(ts);
        (groupType, elementOrder, monogenics) = ComputeDetails(tmpElements);
        elements = new(tmpElements);
    }
    public override T Neutral() => ControlGroup.Neutral();
    public override T Invert(T a) => ControlGroup.Invert(a);
    public override T Op(T a, T b) => ControlGroup.Op(a, b);
    public T Pow(T a, int k)
    {
        if (k == 0)
            return this.Neutral();

        var a0 = k > 0 ? a : this.Invert(a);
        var k0 = k > 0 ? k : -k;
        var acc = a0;
        for (int i = 1; i < k0; ++i)
            acc = this.Op(acc, a0);

        return acc;
    }
    public WorkGroup<T> GenerateProperSubgroup(params T[] ts) => new(this, ts);
    public WorkGroup<T> GenerateSubgroup(IEnumerable<T> ts) => new(this, ts, 0);
    public bool VerifySubGroup(IEnumerable<T> ts)
    {
        if (ts.Any(t => !elements.Contains(t)))
            return false;

        var gens = ComputeGenerators(ts);
        var ts0 = gens.SelectMany(p => p.Value).Select(p => p.e).ToHashSet();
        if (!ts0.SetEquals(ts))
            return false;

        var head = gens.Keys.Select(p => p.e).Ascending().ToHashSet();
        foreach (var e0 in head)
        {
            var ei = Invert(e0);
            if (!ts0.Contains(ei))
                return false;

            foreach (var e1 in head)
            {
                if (!ts0.Contains(this.Op(ei, e1)))
                    return false;
            }
        }
        return true;
    }
}
