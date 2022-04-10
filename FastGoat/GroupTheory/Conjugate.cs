using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{

    public class Normalize<U> : SubGroup<U> where U : struct, IElt
    {
        public Normalize(SubGroup<U> g, SubGroup<U> h, string name, string fmt) : base(g.UpperGroup, name, fmt)
        {
            if (!g.UpperGroup.Equals(h.UpperGroup))
                return;

            Infos = $"in {UpperGroup.Name}";
            Generate(g.AllElements().ToHashSet(new Eq<U>()), h.AllElements().ToHashSet(new Eq<U>()));
        }

        public Normalize(SubGroup<U> g, SubGroup<U> h, string name) : this(g, h, name, g.UpperSet.Fmt) { }
        public Normalize(SubGroup<U> g, SubGroup<U> h) : this(g, h, $"{g.Name}{h.Name}", g.UpperSet.Fmt) { }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

        U Fct(U a, U b) => Op(a, Op(b, Invert(a)));
        void Generate(HashSet<U> g, HashSet<U> h)
        {
            foreach (var e0 in h)
            {
                if (g.All(e1 => h.Contains(Fct(e1, e0))))
                    Add(e0);
            }
        }
    }

    public class Centerize<U> : SubGroup<U> where U : struct, IElt
    {
        public Centerize(SubGroup<U> g, SubGroup<U> h, string name, string fmt) : base(g.UpperGroup, name, fmt)
        {
            if (!g.UpperGroup.Equals(h.UpperGroup))
                return;

            Infos = $"in {UpperGroup.Name}";
            Generate(g.AllElements().ToHashSet(new Eq<U>()), h.AllElements().ToHashSet(new Eq<U>()));
        }

        public Centerize(SubGroup<U> g, SubGroup<U> h, string name) : this(g, h, name, g.UpperSet.Fmt) { }
        public Centerize(SubGroup<U> g, SubGroup<U> h) : this(g, h, $"{g.Name}{h.Name}", g.UpperSet.Fmt) { }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

        U Fct(U a, U b) => Op(a, Op(b, Invert(a)));
        void Generate(HashSet<U> g, HashSet<U> h)
        {
            foreach (var e0 in h)
            {
                if (g.All(e1 => Fct(e1, e0).Equals(e0)))
                    Add(e0);
            }
        }
    }
}
