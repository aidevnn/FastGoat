
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Subgroups;

public struct SubgroupConjugates<T> : IElt<SubgroupConjugates<T>> where T : struct, IElt<T>
{
    public ConcreteGroup<T> Parent { get; }
    public List<ConcreteGroup<T>> Conjugates { get; }
    public ConcreteGroup<T> Representative { get; }
    public int Index { get; }
    public int Order { get; }
    public Dictionary<int, int> Factors { get; }
    public string Subscript { get; set; } = "";
    
    public SubgroupConjugates(ConcreteGroup<T> parent, ConcreteGroup<T> subGroup)
    {
        Parent = parent;
        Conjugates = Group.SubGroupsConjugates(parent, subGroup);
        Representative = Conjugates[0];
        Order = Representative.Count();
        Index = parent.Count() / Order;
        Hash = (Parent.Hash, Order, Index).GetHashCode();
        if (Order == 1)
        {
            Representative.Name = "C1";
            Factors = new() { [1] = 1 };
        }
        else
            Factors = IntExt.PrimesDec(Order);
    }
    
    public SubgroupConjugates(ConcreteGroup<T> parent, List<ConcreteGroup<T>> conjugates)
    {
        Parent = parent;
        Conjugates = conjugates.Count == 0 ? new() { Group.Generate("()", parent, parent.Neutral()) } : conjugates;
        Representative = Conjugates[0];
        Order = Representative.Count();
        Index = parent.Count() / Order;
        Hash = (Parent.Hash, Order, Index).GetHashCode();
        if (Order == 1)
        {
            Representative.Name = "C1";
            Factors = new() { [1] = 1 };
        }
        else
            Factors = IntExt.PrimesDec(Order);
    }
    
    public SubgroupConjugates<T>[] Restriction(ConcreteGroup<T> g)
    {
        if (!g.SubSetOf(Parent))
            throw new GroupException(GroupExceptionType.NotSubGroup);
        
        var lt = Conjugates.Where(h => h.SubSetOf(g)).ToHashSet(new GroupSetEquality<T>());
        var all = new List<SubgroupConjugates<T>>(Size);
        var gens = g.GetGenerators().ToHashSet();
        var act = Group.ByConjugateSet(g);
        while (lt.Count != 0)
        {
            var sg = lt.First();
            var sgSet = sg.ToSet();
            var conjs = Group.Orbits(gens, act, sgSet);
            var subConjs = lt.Where(c0 => conjs.Any(c1 => c0.SetEquals(c1))).ToList();
            all.Add(new SubgroupConjugates<T>(g, subConjs));
            lt.ExceptWith(subConjs);
        }
        
        return all.ToArray();
    }
    
    public bool IsPGroup()
    {
        return Factors.Count == 1;
    }

    public bool Contains(ConcreteGroup<T> g) => Conjugates.Any(e => e.SetEquals(g));
    public bool Contains(HashSet<T> g) => Conjugates.Any(e => e.SetEquals(g));
    public bool SubSetOf(ConcreteGroup<T> g) => Conjugates.All(e => g.SuperSetOf(e));
    public GroupType GroupType => Representative.GroupType;
    public bool IsMonogenic => Representative.GetGenerators().Count() == 1;
    public int Size => Conjugates.Count;
    public (int, int, GroupType) OST => (Order, Size, GroupType);
    public bool IsNormal => Size == 1;
    public bool IsProper => Index != 1;
    public bool IsTrivial => Order == 1;
    public bool IsProperNormal => IsNormal && IsProper;

    public bool IsSubClassOf(SubgroupConjugates<T> other) => Order < other.Order && 
                                                             Conjugates.Any(cj => cj.SubSetOf(other.Representative));
    public bool IsSuperClassOf(SubgroupConjugates<T> other) => Order > other.Order && 
                                                               Conjugates.Any(cj => cj.SuperSetOf(other.Representative));
    public bool Equals(SubgroupConjugates<T> other) => Hash == other.Hash && 
                                                       Parent.SetEquals(other.Parent) &&
                                                       Conjugates.Any(e => e.SetEquals(other.Representative));
    
    public int CompareTo(SubgroupConjugates<T> other) => OST.CompareTo(other.OST);

    public int Hash { get; }

    public SubgroupConjugates<WElt> ToGroupWrapper()
    {
        if (this is SubgroupConjugates<WElt> cj)
            return cj;
        
        var gt = Parent.ToGroupWrapper();
        var sub = Group.Generate(Representative.Name, gt, Representative.GetGenerators().Select(e => new WElt(e)).ToArray());
        return new(gt, sub);
    }

    public string FullName => Subscript.Length == 0 ? Representative.Name : $"{Representative.NameParenthesis()}{Subscript}";

    public override int GetHashCode() => Hash;
    public override string ToString() => Representative.Name;
}