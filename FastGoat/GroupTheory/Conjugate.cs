using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{

    public class Normalize<U> : SubGroup<U> where U : struct, IElt
    {
        public Normalize(SubGroup<U> h, SubSet<U> s) : base(h.UpperGroup)
        {
            if (!h.UpperGroup.Equals(s.UpperSet))
                return;

            if (!h.IsGroup)
                return;

            Infos = $"in {UpperGroup.Name}";
            Generate(h.AllElements().ToHashSet(new Eq<U>()), s.AllElements().ToHashSet(new Eq<U>()));
            Name = "K";
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

        void Generate(HashSet<U> h, HashSet<U> s)
        {
            foreach (var x in h)
            {
                if (s.All(e1 => s.Contains(Op(Op(x, e1), Invert(x)))))
                    Add(x);
            }
        }
    }

    public class Centerize<U> : SubGroup<U> where U : struct, IElt
    {
        public Centerize(SubGroup<U> h, SubSet<U> s) : base(h.UpperGroup)
        {
            if (!h.UpperGroup.Equals(s.UpperSet))
                return;

            Infos = $"in {UpperGroup.Name}";
            Generate(h.AllElements().ToHashSet(new Eq<U>()), s.AllElements().ToHashSet(new Eq<U>()));
            Name = "Z";
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a) => UpperGroup.Invert(a);
        public override U Op(U a, U b) => UpperGroup.Op(a, b);

        void Generate(HashSet<U> h, HashSet<U> s)
        {
            foreach (var x in h)
            {
                if (s.All(e1 => e1.Equals(Op(Op(x, e1), Invert(x)))))
                    Add(x);
            }
        }
    }
}
