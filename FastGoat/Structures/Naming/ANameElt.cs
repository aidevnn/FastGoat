using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Naming;

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
    public (GroupType, int, int, int, int, string) Infos => (ContentGroup!.GroupType, NbExt, NbSdp, NbDp, Name.Length, Name);
    public int Hash => (ContentType, Name).GetHashCode();
    public string NameParenthesis => Name.WithParenthesis();
    public bool Equals(ANameElt? other) => other?.Name == Name;

    public int CompareTo(ANameElt? other)
    {
        if (other is null)
            return 1;

        return Infos.CompareTo(other.Infos);
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}