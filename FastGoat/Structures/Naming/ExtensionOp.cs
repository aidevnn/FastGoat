using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Naming;

public class ExtensionOp : ANameElt
{
    public ANameElt Lhs { get; }
    public ANameElt Rhs { get; }

    public ExtensionOp(ANameElt lhs, ANameElt rhs, ConcreteGroup<TableElt> g)
    {
        ContentGroup = g;
        ContentType = NodeType.Extension;
        (Lhs, Rhs) = (lhs, rhs);
        Name = $"{Lhs.NameParenthesis} . {Rhs.NameParenthesis}";
        Depth = 1 + int.Max(Lhs.Depth, Rhs.Depth);
        Weight = Lhs.Weight * Rhs.Weight;
    }
}