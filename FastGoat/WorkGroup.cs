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
    private WorkGroup(WorkGroup<T> group, IEnumerable<T> ts, bool verify) : base(group.BaseGroup)
    {
        ControlGroup = group;
        if (ts.Any(e => !group.BaseGroup.Equals(e.Group)))
            throw new BaseGroupException();

        if (!verify)
        {
            var tmpElements = Generate(ts);
            (groupType, elementOrder, monogenics) = ComputeDetails(tmpElements);
            elements = new(tmpElements);
        }
        else
        {
            var (sgENum, gens) = VerifySubGroup(ts);
            if (sgENum == SubGroupEnum.IsSubGroup)
            {
                var tmpElements = ts;
                (groupType, elementOrder, monogenics) = ComputeDetails(gens);
                elements = new(tmpElements);
            }
            else if (sgENum == SubGroupEnum.OutOfSubGroup)
                throw new SubGroupException("Set is not a valid subGroup");
            else
                throw new BaseGroupException();
        }
    }
    public override T Neutral() => ControlGroup.Neutral();
    public override T Invert(T a) => ControlGroup.Invert(a);
    public override T Op(T a, T b) => ControlGroup.Op(a, b);
    public WorkGroup<T> GenerateProperSubgroup(params T[] ts) => new(this, ts);
    public WorkGroup<T> GenerateSubgroup(IEnumerable<T> ts) => new(this, ts, false);
    public WorkGroup<T> GetSubgroup(IEnumerable<T> ts) => new(this, ts, true);
}
