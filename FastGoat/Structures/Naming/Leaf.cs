using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Naming;

public class Leaf : ANameElt
{
    public Leaf(AllSubgroups<TableElt> subgroups, DecompType decompType)
    {
        ContentGroup = subgroups.Parent;
        ContentType = NodeType.Leaf;
        Depth = 0;
        if (decompType == DecompType.Abelian)
        {
            var abType = Group.AbelianGroupType(ContentGroup);
            Name = abType.Glue(" x ", "C{0}");
            Weight = abType.Length * 5 + 15;
        }
        else if (decompType == DecompType.SimpleNonAbelian)
        {
            Name = SimpleNonAbelians(ContentGroup);
            Weight = 30;
        }
    }

    public Leaf(ConcreteGroup<TableElt> g)
    {
        ContentGroup = g;
        if (g.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.GroupDef);

        ContentType = NodeType.Leaf;
        var abType = Group.AbelianGroupType(ContentGroup);
        Name = abType.Glue(" x ", "C{0}");
        Weight = abType.Length * 5 + 15;
        Depth = 0;
    }

    public Leaf(ConcreteGroup<TableElt> g, string name)
    {
        ContentGroup = g;
        ContentType = NodeType.Leaf;
        Name = name;
        Depth = 0;
        Weight = 30;
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