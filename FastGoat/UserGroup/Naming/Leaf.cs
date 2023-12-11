using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.Naming;

public class Leaf : ANameElt
{
    public Leaf(AllSubgroups<TableElt> subgroups, DecompType decompType)
    {
        ContentGroup = subgroups.Parent;
        ContentType = NodeType.Leaf;
        if (decompType == DecompType.Abelian)
            Name = Group.AbelianInvariants(ContentGroup).Select(e => e.o).Glue(" x ", "C{0}");
        else if (decompType == DecompType.SimpleNonAbelian)
            Name = SimpleNonAbelians(ContentGroup);
    }

    public Leaf(ConcreteGroup<TableElt> g)
    {
        ContentGroup = g;
        if (g.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.GroupDef);

        ContentType = NodeType.Leaf;
        Name = Group.AbelianInvariants(ContentGroup).Select(e => e.o).Glue(" x ", "C{0}");
    }

    public Leaf(ConcreteGroup<TableElt> g, string name)
    {
        ContentGroup = g;
        ContentType = NodeType.Leaf;
        Name = name;
    }
    
    static string SimpleNonAbelians(ConcreteGroup<TableElt> G)
    {
        var og = G.Count();
        if (og == 60)
            return "A5";
        if (og == 168)
            return "SL(3,2)";
        if (og == 360)
            return "A6";
        if (og == 504)
            return "SL(2,8)";
        if (og == 660)
            return "L2(11)";
        if (og == 1092)
            return "L2(13)";

        return G.Name;
    }
}