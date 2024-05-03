using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Naming;

public class DirectProductOp : ANameElt
{
    public ANameElt[] Elts { get; }

    public DirectProductOp(ANameElt lhs, ANameElt rhs, ConcreteGroup<WElt> g)
    {
        ContentType = NodeType.DirectProduct;
        ContentGroup = g;
        if (lhs.ContentType == NodeType.DirectProduct && rhs.ContentType == NodeType.DirectProduct)
            Elts = ((DirectProductOp)lhs).Elts.Concat(((DirectProductOp)rhs).Elts).ToArray();
        else if (lhs.ContentType == NodeType.DirectProduct)
            Elts = ((DirectProductOp)lhs).Elts.Append(rhs).ToArray();
        else if (rhs.ContentType == NodeType.DirectProduct)
            Elts = ((DirectProductOp)rhs).Elts.Append(lhs).ToArray();
        else
            Elts = new[] {lhs, rhs};

        var abGens = Elts.Where(e => e.ContentGroup!.GroupType == GroupType.AbelianGroup)
            .SelectMany(e => e.ContentGroup!.GetGenerators()).Distinct().ToArray();
        if (abGens.Length != 0)
        {
            var nab = Elts.Where(e => e.ContentGroup!.GroupType == GroupType.NonAbelianGroup).ToArray();
            var ab = Group.Generate("Ab", g, abGens);
            var leaf = new Leaf(ab);
            Elts = nab.Prepend(leaf).ToArray();
        }

        Elts = Elts.Order().ToArray();
        Name = Elts.Select(e => e.NameParenthesis).Glue(" x ");
        Depth = 1 + Elts.Max(e => e.Depth);
        Weight = Elts.Max(e => e.Weight);
    }
}