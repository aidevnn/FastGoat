using System;
using System.Collections.Generic;
using System.Linq;

namespace FastGoat.SetTheory
{
    public class EmptySet<U> : SubSet<U> where U : struct, IElt
    {
        public EmptySet(FSet<U> fset, string name, string fmt) : base(fset, name, fmt) { }
        public EmptySet(FSet<U> fset, string name) : base(fset, name, fset.Fmt) { }
        public EmptySet(FSet<U> fset) : base(fset, "{}", fset.Fmt) { }
    }

    public class Inter<U> : SubSet<U> where U : struct, IElt
    {
        public Inter(SubSet<U> g, SubSet<U> h, string name, string fmt) : base(g.UpperSet, name, fmt)
        {
            if (!g.UpperSet.Equals(h.UpperSet))
                return;

            foreach (var e in h.AllElements().Intersect(g.AllElements()))
                Add(e);

            Infos = $"in {UpperSet.Name}";
        }

        public Inter(SubSet<U> g, SubSet<U> h, string name) : this(g, h, name, g.UpperSet.Fmt) { }
        public Inter(SubSet<U> g, SubSet<U> h) : this(g, h, $"{g.Name}_{h.Name}", g.UpperSet.Fmt) { }
    }

    public class Union<U> : SubSet<U> where U : struct, IElt
    {
        public Union(SubSet<U> g, SubSet<U> h, string name, string fmt) : base(g.UpperSet, name, fmt)
        {
            if (!g.UpperSet.Equals(h.UpperSet))
                return;

            foreach (var e in h.AllElements().Union(g.AllElements()))
                Add(e);

            Infos = $"in {UpperSet.Name}";
        }

        public Union(SubSet<U> g, SubSet<U> h, string name) : this(g, h, name, g.UpperSet.Fmt) { }
        public Union(SubSet<U> g, SubSet<U> h) : this(g, h, $"{g.Name}_{h.Name}", g.UpperSet.Fmt) { }

        public Union(FSet<U> fset, string name, string fmt, params U[] us) : base(fset, name, fmt)
        {
            if (us.Any(e => !UpperSet.Equals(e.FSet)))
                return;

            foreach (var e in us)
                Add(e);

            Infos = $"in {UpperSet.Name}";
        }

        public Union(FSet<U> fset, string name, params U[] us) : this(fset, name, fset.Fmt, us) { }
        public Union(FSet<U> fset, params U[] us) : this(fset, "G", fset.Fmt, us) { }
    }

}
