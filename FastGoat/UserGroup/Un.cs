using System.CodeDom;

namespace FastGoat.UserGroup;

public class Un : ConcreteGroup<Automorphism<ZnInt>>
{
    public Un(int n) : base($"U{n}", Group.Aut(new Cn(n)), true)
    {
        var autCn = (AutomorphismGroup<ZnInt>)BaseGroup;
        var cn = autCn.G;
        var elements = new List<Automorphism<ZnInt>>();
        for (int k = 1; k < n; ++k)
        {
            if (IntExt.Gcd(k, n) != 1)
                continue;

            var ak = autCn[(cn[1], cn[k])];
            elements.Add(ak);
        }

        Hash = (BaseGroup.Hash, "Un").GetHashCode();
        Elements = elements.ToHashSet();
        LongestCycles = Group.LongestCycles(autCn, elements);
        ElementsOrders = Group.ElementsOrders(LongestCycles);
        GroupType = Group.IsCommutative(autCn, LongestCycles.Keys) ? GroupType.AbelianGroup : GroupType.NonAbelianGroup;
    }
}