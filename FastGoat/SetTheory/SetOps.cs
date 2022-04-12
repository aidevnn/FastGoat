using System;
using System.Collections.Generic;
using System.Linq;

namespace FastGoat.SetTheory
{
    public class EmptySet<U> : SubSet<U> where U : struct, IElt
    {
        public EmptySet(FSet<U> fset) : base(fset) { Name = "{}"; }
    }

    public class Inter<U> : SubSet<U> where U : struct, IElt
    {
        public Inter(SubSet<U> g, SubSet<U> h) : base(g.UpperSet)
        {
            if (!g.UpperSet.Equals(h.UpperSet))
                return;

            foreach (var e in h.AllElements().Intersect(g.AllElements()))
                Add(e);

            Infos = $"in {UpperSet.Name}";
            Name = $"{g.Name}_{h.Name}";
        }
    }

    public class Union<U> : SubSet<U> where U : struct, IElt
    {
        public Union(SubSet<U> g, SubSet<U> h) : base(g.UpperSet)
        {
            if (!g.UpperSet.Equals(h.UpperSet))
                return;

            foreach (var e in h.AllElements().Union(g.AllElements()))
                Add(e);

            Infos = $"in {UpperSet.Name}";
            Name = $"{g.Name}_{h.Name}";
        }

        public Union(FSet<U> fset, params U[] us) : base(fset)
        {
            if (us.Any(e => !UpperSet.Equals(e.FSet)))
                return;

            foreach (var e in us)
                Add(e);

            Infos = $"in {UpperSet.Name}";
            Name = "G";
        }
    }

}
