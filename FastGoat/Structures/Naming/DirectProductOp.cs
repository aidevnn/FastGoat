using System.Numerics;
using FastGoat.Commons;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.Structures.Naming;

public class DirectProductOp : ANameElt
{
    public ANameElt[] Elts { get; }

    public DirectProductOp(ANameElt lhs, ANameElt rhs, ConcreteGroup<TableElt> g)
    {
        ContentType = NodeType.DirectProduct;
        ContentGroup = g;
        if (lhs.ContentType == NodeType.DirectProduct && rhs.ContentType == NodeType.DirectProduct)
            Elts = [..((DirectProductOp)lhs).Elts, ..((DirectProductOp)rhs).Elts];
        else if (lhs.ContentType == NodeType.DirectProduct)
            Elts = [..((DirectProductOp)lhs).Elts, rhs];
        else if (rhs.ContentType == NodeType.DirectProduct)
            Elts = [..((DirectProductOp)rhs).Elts, lhs];
        else
            Elts = [lhs, rhs];

        var abGens = Elts.Where(e => e.ContentGroup!.GroupType == GroupType.AbelianGroup)
            .SelectMany(e => e.ContentGroup!.GetGenerators()).Distinct().ToArray();
        if (abGens.Length != 0)
        {
            var nab = Elts.Where(e => e.ContentGroup!.GroupType == GroupType.NonAbelianGroup).ToArray();
            var ab = Group.Generate("Ab", g, abGens);
            var leaf = new Leaf(ab);
            Elts = [leaf, ..nab];
        }

        Elts = Elts.Order().ToArray();
        Name = Elts.Select(e => e.NameParenthesis).Glue(" x ");
        Depth = 1 + Elts.Max(e => e.Depth);
        Weight = Elts.Aggregate(BigInteger.One, (a, b) => a * b.Weight);
    }
}