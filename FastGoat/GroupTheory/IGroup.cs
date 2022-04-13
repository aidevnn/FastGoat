using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public interface IGroup<U> : ISubSet<U> where U : struct, IElt
    {
        U Neutral { get; }
        U Invert(U a);
        U Op(U a, U b);
    }

    public abstract class Group<U> : FSet<U>, IGroup<U> where U : struct, IElt
    {
        protected Group(int capacity) : base(capacity) { }

        public FSet<U> UpperSet => this;

        public abstract U Neutral { get; }
        public abstract U Invert(U a);
        public abstract U Op(U a, U b);

        public bool IsEqual(SubSet<U> set)
        {
            if (!UpperSet.Equals(set.UpperSet)) return false;
            return Elts.SetEquals(set.AllElements());
        }
    }

    public interface ISubGroup<U> : IGroup<U> where U : struct, IElt
    {
        Group<U> UpperGroup { get; }
    }
}
