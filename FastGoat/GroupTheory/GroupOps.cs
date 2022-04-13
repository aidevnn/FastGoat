using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public class Monogenic<U> : SubGroup<U> where U : struct, IElt
    {
        public Monogenic(Group<U> group, U e) : base(group)
        {
            Infos = $"in {group.Name}";
            Name = $"<{e}>";
            Add(group.Neutral);
            if (group.Equals(e.FSet))
                Generate(e);
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

        void Generate(U e)
        {
            var acc = e;
            while (!acc.Equals(Neutral))
            {
                Add(acc);
                acc = Op(e, acc);
            }
        }
    }

    public class DirectProduct<U> : SubGroup<U> where U : struct, IElt
    {
        public DirectProduct(SubGroup<U> g, SubGroup<U> h) : base(g.UpperGroup)
        {
            if (!g.UpperGroup.Equals(h.UpperGroup))
                return;

            foreach (var e in h.AllElements().Union(g.AllElements()))
                Add(e);

            Infos = $"in {UpperGroup.Name}";
            Name = $"{g.Name}.{h.Name}";
            Generate();
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

        void Generate()
        {
            var elts = Elts.ToList();
            foreach (var e0 in elts)
            {
                foreach (var e1 in elts)
                {
                    var e2 = Op(e0, e1);
                    Add(e2);
                }
            }
        }
    }

    public class SubGroupOf<U> : SubGroup<U> where U : struct, IElt
    {
        public SubGroupOf(Group<U> group, SubSet<U> subSet) :base(group) 
        {
            if (!group.Equals(subSet.UpperSet))
                return;

            if (subSet.AllElements().Any(e => !group.Contains(e)))
                return;

            foreach (var e in subSet.AllElements())
                Add(e);
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

    }

    public class GroupOp<U> : SubGroup<U> where U : struct, IElt
    {
        public GroupOp(SubGroup<U> sub, U e) : base(sub.UpperGroup)
        {
            if (!sub.UpperSet.Equals(e.FSet))
                return;

            Generate(sub.AllElements(), e);
            Name = "Hx";
        }

        public GroupOp(U e, SubGroup<U> sub) : base(sub.UpperGroup)
        {
            if (!sub.UpperSet.Equals(e.FSet))
                return;

            Generate(e, sub.AllElements());
            Name = "xH";
        }

        void Generate(U e, IEnumerable<U> sub)
        {
            foreach (var e0 in sub)
                Add(Op(e, e0));
        }

        void Generate(IEnumerable<U> sub, U e)
        {
            foreach (var e0 in sub)
                Add(Op(e0, e));
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);
    }
}
