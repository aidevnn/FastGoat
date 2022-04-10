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
            return new SubGroupOf<U>(group, subSet, subSet.Name);
        }

        public static SubGroup<U> SubGroupOf<U>(this SubSet<U> subSet, Group<U> group, string name) where U : struct, IElt
        {
            return new SubGroupOf<U>(group, subSet, name);
        }

        public static SubGroup<U> Monogenic<U>(this Group<U> group, U e) where U : struct, IElt
        {
            return new Monogenic<U>(group, e);
        }

        public static SubGroup<U> Monogenic<U>(this Group<U> group, U e, string name) where U : struct, IElt
        {
            return new Monogenic<U>(group, e, name);
        }

        public static SubGroup<U> DirectProduct<U>(this SubGroup<U> g, SubGroup<U> h) where U : struct, IElt
        {
            return new DirectProduct<U>(g, h);
        }

        public static SubGroup<U> DirectProduct<U>(this SubGroup<U> g, SubGroup<U> h, string name) where U : struct, IElt
        {
            return new DirectProduct<U>(g, h, name);
        }

        public static SubGroup<U> Union<U>(this Group<U> g, params U[] us) where U : struct, IElt
        {
            return g.EmptySet().Union(us).SubGroupOf(g);
        }

        public static SubGroup<U> Union<U>(this Group<U> g, string name, params U[] us) where U : struct, IElt
        {
            return g.EmptySet().Union(us).SubGroupOf(g, name);
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

        public static SubGroup<U> Develop<U>(this SubGroup<U> subGroup, string name) where U : struct, IElt
        {
            var g = subGroup.Develop();
            g.Name = name;
            return g;
        }

        public static SubGroup<U> Develop<U>(this SubSet<U> subSet, Group<U> group) where U : struct, IElt
        {
            return subSet.SubGroupOf(group).Develop();
        }
        public static SubGroup<U> Develop<U>(this SubSet<U> subSet, Group<U> group, string name) where U : struct, IElt
        {
            return subSet.SubGroupOf(group).Develop(name);
        }

        public static QuotientGroup<U> Over<U>(this SubGroup<U> g, SubGroup<U> h) where U : struct, IElt
        {
            return new QuotientGroup<U>(g, h);
        }
    }
}
