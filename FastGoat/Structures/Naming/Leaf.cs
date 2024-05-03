using System.Numerics;
using System.Text.RegularExpressions;
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Subgroups;

namespace FastGoat.Structures.Naming;

public class Leaf : ANameElt
{
    public Leaf(AllSubgroups<WElt> subgroups, DecompType decompType)
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

    public Leaf(ConcreteGroup<WElt> g)
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

    public Leaf(ConcreteGroup<WElt> g, string name)
    {
        ContentGroup = g;
        ContentType = NodeType.Leaf;
        Name = name;
        Depth = 0;
        Weight = 30;
    }

    static string SimpleNonAbelians(ConcreteGroup<WElt> G)
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
        if (og == 2520)
            return "A7";

        return G.Name;
    }

    public (string prefix, int[] coefs) LeafDetails()
    {
        if (ContentGroup!.GroupType == GroupType.AbelianGroup)
            return ("C", Group.AbelianGroupType(ContentGroup));

        string regX = @"(Q|Dic|A|S|L2)(\d+)|(GL|SL)\((\d+),(\d+)\)";
        var match = Regex.Match(Name, regX);
        var s1 = match.Groups["1"].Value;
        var s3 = match.Groups["3"].Value;
        if (string.IsNullOrEmpty(s3) && !string.IsNullOrEmpty(s1))
        {
            var n = int.Parse(match.Groups["2"].Value);
            return (s1, new[] { n });
        }
        else if (string.IsNullOrEmpty(s1) && !string.IsNullOrEmpty(s3))
        {
            var n = int.Parse(match.Groups["4"].Value);
            var p = int.Parse(match.Groups["5"].Value);
            return (s3, new[] { n, p });
        }

        throw new();
    }
}