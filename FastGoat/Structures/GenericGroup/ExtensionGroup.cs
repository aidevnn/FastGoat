using FastGoat.Structures.CartesianProduct;

namespace FastGoat.Structures.GenericGroup;

public class ExtensionGroup<Tn, Tg> : ConcreteGroup<Ep2<Tn, Tg>> where Tg : struct, IElt<Tg> where Tn : struct, IElt<Tn>
{
    public MapElt<Ep2<Tg, Tg>, Tn> Map => ExtBase.Map;
    public MapElt<Tg, Automorphism<Tn>> L => ExtBase.L;
    public ConcreteGroup<Tg> G => ExtBase.G;
    public ConcreteGroup<Tn> N => ExtBase.N;

    public ExtensionGroupBase<Tn, Tg> ExtBase { get; }

    public ExtensionGroup(ConcreteGroup<Tn> n, MapElt<Tg, Automorphism<Tn>> l, MapElt<Ep2<Tg, Tg>, Tn> map, ConcreteGroup<Tg> g)
        : this(new ExtensionGroupBase<Tn, Tg>(n, l, map, g))
    {
    }

    private ExtensionGroup(ExtensionGroupBase<Tn, Tg> ExtBase) : base("Ext", ExtBase)
    {
        this.ExtBase = ExtBase;
        var nName = ExtBase.N.Name.Contains('x') ? $"({N})" : $"{N}";
        var gName = ExtBase.G.Name.Contains('x') ? $"({G})" : $"{G}";
        Name = $"{nName} . {gName}";
        Hash = (Name, Elements.Count).GetHashCode();
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
}