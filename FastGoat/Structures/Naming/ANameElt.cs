using System.Numerics;
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

    public BigInteger Weight { get; set; }
    public int Depth { get; set; }

    public ConcreteGroup<TableElt>? ContentGroup { get; set; }
    public NodeType ContentType { get; set; }
    public string Name { get; set; } = "C1";
    public (int, string) LN => (Name.Length, Name);
    public int Hash => (ContentType, Name).GetHashCode();
    public string NameParenthesis => Name.WithParenthesis();
    public bool Equals(ANameElt? other) => other?.Name == Name;

    public bool IsExtension => ContentType == NodeType.Extension;

    public int CompareTo(ANameElt? other)
    {
        if (other is null)
            return 1;

        var compCTC = IsExtension.CompareTo(other.IsExtension);
        if (compCTC != 0)
            return compCTC;

        var compD = Depth.CompareTo(other.Depth);
        if (compD != 0)
            return compD;

        var compCT = ContentType.CompareTo(other.ContentType);
        if (compCT != 0)
            return compCT;

        var compW = Weight.CompareTo(other.Weight);
        if (compW != 0)
            return compW;

        return LN.CompareTo(other.LN);
    }

    public override int GetHashCode() => Hash;

    public override string ToString() => Name;
}