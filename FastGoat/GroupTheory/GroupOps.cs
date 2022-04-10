using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public class Monogenic<U> : SubGroup<U> where U : struct, IElt
    {
        public Monogenic(Group<U> group, U e, string name, string fmt) : base(group, name, fmt)
        {
            Infos = $"in {group.Name}";
            Add(group.Neutral);
            if (group.Equals(e.FSet))
                Generate(e);
        }

        public Monogenic(Group<U> group, U e, string name) : this(group, e, name, group.Fmt) { }
        public Monogenic(Group<U> group, U e) : this(group, e, $"<{e}>", group.Fmt) { }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

        void Generate(U e)
        {
            var acc = e;
            while (!acc.Equals(Neutral))
            {
                Add(acc);
                if (OpLR == XOpLR.Left)
                    acc = Op(e, acc);
                else
                    acc = Op(acc, e);
            }
        }
    }

    public class DirectProduct<U> : SubGroup<U> where U : struct, IElt
    {
        public DirectProduct(SubGroup<U> g, SubGroup<U> h, string name, string fmt) : base(g.UpperGroup, name, fmt)
        {
            if (!g.UpperGroup.Equals(h.UpperGroup))
                return;

            foreach (var e in h.AllElements().Union(g.AllElements()))
                Add(e);

            Infos = $"in {UpperGroup.Name}";
            Generate();
        }

        public DirectProduct(SubGroup<U> g, SubGroup<U> h, string name) : this(g, h, name, g.UpperSet.Fmt) { }
        public DirectProduct(SubGroup<U> g, SubGroup<U> h) : this(g, h, $"{g.Name}.{h.Name}", g.UpperSet.Fmt) { }

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
        public SubGroupOf(Group<U> group, SubSet<U> subSet, string name, string fmt) :base(group,name,fmt) 
        {
            if (!group.Equals(subSet.UpperSet))
                return;

            if (subSet.AllElements().Any(e => !group.Contains(e)))
                return;

            foreach (var e in subSet.AllElements())
                Add(e);
        }

        public SubGroupOf(Group<U> group, SubSet<U> subSet, string name) : this(group, subSet, name, group.Fmt) { }
        public SubGroupOf(Group<U> group, SubSet<U> subSet) : this(group, subSet, "G", group.Fmt) { }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

    }

    public class GroupOp<U> : SubGroup<U> where U : struct, IElt
    {
        public GroupOp(SubGroup<U> sub, U e, string name, string fmt) : base(sub.UpperGroup, name, fmt)
        {
            OpLR = XOpLR.Right;
            if (!sub.UpperSet.Equals(e.FSet))
                return;

            Generate(sub.AllElements().ToArray(), e);
        }

        public GroupOp(SubGroup<U> sub, U e, string name) : this(sub, e, name, sub.UpperSet.Fmt) { }
        public GroupOp(SubGroup<U> sub, U e) : this(sub, e, "Hx", sub.UpperSet.Fmt) { }

        public GroupOp(U e, SubGroup<U> sub, string name, string fmt) : base(sub.UpperGroup, name, fmt)
        {
            OpLR = XOpLR.Left;
            if (!sub.UpperSet.Equals(e.FSet))
                return;

            Generate(sub.AllElements().ToArray(), e);
        }

        public GroupOp(U e, SubGroup<U> sub, string name) : this(e, sub, name, sub.UpperSet.Fmt) { }
        public GroupOp(U e, SubGroup<U> sub) : this(e, sub, "xH", sub.UpperSet.Fmt) { }

        void Generate(U[] sub, U e)
        {
            if (OpLR == XOpLR.Right)
            {
                foreach (var e0 in sub)
                    Add(Op(e, e0));
            }
            else
            {
                foreach (var e0 in sub)
                    Add(Op(e0, e));
            }
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

    }

}
