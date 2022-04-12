using System;
using System.Collections.Generic;
using System.Linq;

namespace FastGoat.SetTheory
{
    public static class SetExt
    {
        public static SubSet<U> EmptySet<U>(this FSet<U> fset) where U : struct, IElt
        {
            return new EmptySet<U>(fset);
        }

        public static SubSet<U> Union<U>(this SubSet<U> g, SubSet<U> h) where U : struct, IElt
        {
            return new Union<U>(g, h);
        }

        public static SubSet<U> Union<U>(this SubSet<U> g, params U[] us) where U : struct, IElt
        {
            return new Union<U>(g.UpperSet, g.AllElements().Union(us, new Eq<U>()).ToArray());
        }

        public static SubSet<U> Union<U>(this FSet<U> fset, params U[] us) where U : struct, IElt
        {
            return fset.EmptySet().Union(us);
        }

        public static SubSet<U> Inter<U>(this SubSet<U> g, SubSet<U> h) where U : struct, IElt
        {
            return new Inter<U>(g, h);
        }

    }
}
