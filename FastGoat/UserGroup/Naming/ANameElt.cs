using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.Naming;

public abstract class ANameElt : IElt<ANameElt>
{
    [Flags]
    public enum DecompType
    {
        Abelian = 0,
        DirectProduct = 1,
        SemiDirectProduct = 2,
        Extension = 3,
        SimpleNonAbelian = 4
    }
    
    [Flags]
    public enum NodeType
    {
        Leaf = 0,
        DirectProduct = 1,
        SemiDirectProduct = 2,
        Extension = 3
    }
    
    public ConcreteGroup<TableElt>? ContentGroup { get; set; }
    public NodeType ContentType { get; set; }
    public string Name { get; set; } = "C1";
    public int NbDp => System.Text.RegularExpressions.Regex.Count(Name, @"x ");
    public int NbSdp => System.Text.RegularExpressions.Regex.Count(Name, @"x\:");
    public int NbExt => System.Text.RegularExpressions.Regex.Count(Name, @"\. ");
    public (GroupType, int, int, int, int) Infos => (ContentGroup!.GroupType, NbExt, NbSdp, NbDp, Name.Length);
    public int Hash => (ContentType, Name).GetHashCode();
    public string NameParenthesis => Name.WithParenthesis();
    public bool Equals(ANameElt? other) => other?.Name == Name;

    public int CompareTo(ANameElt? other)
    {
        if (other is null)
            return 1;

        var comp = Infos.CompareTo(other.Infos);
        if (comp != 0)
            return comp;

        return string.Compare(Name, other.Name, StringComparison.CurrentCulture);
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}