using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;

namespace FastGoat.UserGroup.Integers;

public class Un : ConcreteGroup<Automorphism<ZnInt>>
{
    public Cn Cn { get; }

    public Un(int n) : base($"U{n}", Group.AutBase(new Cn(n)), true)
    {
        var autCn = (AutomorphismGroup<ZnInt>)BaseGroup;
        Cn = (Cn)autCn.G;
        var elements = new List<Automorphism<ZnInt>>();
        for (int k = 1; k < n; ++k)
        {
            if (IntExt.Gcd(k, n) != 1)
                continue;

            var ak = autCn[(Cn[1], Cn[k])];
            elements.Add(ak);
        }

        Hash = (BaseGroup.Hash, "Un").GetHashCode();
        Elements = elements.ToHashSet();
        ElementsOrders = Group.ElementsOrders(autCn, Elements);
        var (tmpElements, uniqueGenerators) = Group.UniqueGenerators(this, Elements.ToArray());
        PseudoGenerators = new(uniqueGenerators);
        GroupType = Group.IsCommutative(autCn, PseudoGenerators) ? GroupType.AbelianGroup : GroupType.NonAbelianGroup;
    }
}