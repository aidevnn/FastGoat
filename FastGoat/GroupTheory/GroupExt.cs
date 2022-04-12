using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public static class GroupExt
    {
        public static SubGroup<U> SubGroupOf<U>(this SubSet<U> subSet, Group<U> group) where U : struct, IElt
        {
            return new SubGroupOf<U>(group, subSet).GName(subSet.Name);
        }

        public static SubGroup<U> Monogenic<U>(this Group<U> group, U e) where U : struct, IElt
        {
            return new Monogenic<U>(group, e);
        }

        public static SubGroup<U> DirectProduct<U>(this SubGroup<U> g, SubGroup<U> h) where U : struct, IElt
        {
            return new DirectProduct<U>(g, h);
        }

        public static SubGroup<U> Union<U>(this Group<U> g, params U[] us) where U : struct, IElt
        {
            return g.EmptySet().Union(us).SubGroupOf(g);
        }

        public static SubGroup<U> Develop<U>(this SubGroup<U> subGroup) where U : struct, IElt
        {
            SubGroup<U> g = subGroup;
            int sz = 0;
            do
            {
                sz = g.Count;
                g = g.DirectProduct(subGroup);
            } while (sz != g.Count);

            return g;
        }

        public static SubGroup<U> Develop<U>(this SubSet<U> subSet, Group<U> group) where U : struct, IElt
        {
            return subSet.SubGroupOf(group).Develop();
        }

        public static SubGroup<U> GeneratedSubGroup<U>(this Group<U> g, params U[] us) where U : struct, IElt
        {
            return g.EmptySet().Union(us).SubGroupOf(g).Develop();
        }

        public static QuotientGroup<U> Over<U>(this SubGroup<U> g, SubGroup<U> h) where U : struct, IElt
        {
            return new QuotientGroup<U>(g, h);
        }

        public static SubGroup<U> Normalize<U>(this SubGroup<U> g, SubSet<U> h) where U : struct, IElt
        {
            return new Normalize<U>(g, h);
        }
        public static SubGroup<U> Normalize<U>(this SubGroup<U> g, params U[] us) where U : struct, IElt
        {
            return new Normalize<U>(g, g.UpperSet.Union(us));
        }

        public static SubGroup<U> Centerize<U>(this SubGroup<U> g, SubSet<U> h) where U : struct, IElt
        {
            return new Centerize<U>(g, h);
        }
        public static SubGroup<U> Centerize<U>(this SubGroup<U> g, params U[] us) where U : struct, IElt
        {
            return new Centerize<U>(g, g.UpperSet.Union(us));
        }

    }
}
