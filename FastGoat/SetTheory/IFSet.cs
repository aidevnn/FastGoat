using System;
using System.Collections.Generic;
using System.Linq;

namespace FastGoat.SetTheory
{
    public interface IFSet : IEquatable<IFSet>
    {
        string Name { get; set; }
        string Fmt { get; set; }
        string FmtElt { get; set; }
    }

    public interface IElt : IEquatable<IElt>, IComparable<IElt>
    {
        int HashCode { get; }
        int[] Table { get; }
        IFSet FSet { get; }
        string Display(string name = "");
        string EltStr();
        string Infos();
    }

    public interface IFSet<U> : IFSet where U : struct, IElt
    {
        void Add(U e);
        bool Contains(U e);
        int Count { get; }
        IEnumerable<U> AllElements();
    }

    public abstract class FSet<U> : IFSet<U> where U : struct, IElt
    {
        protected FSet(int capacity)
        {
            Capacity = capacity > 1000000 ? 1000000 : capacity;
            Elts = new HashSet<U>(Capacity, new Eq<U>());
        }

        public int Capacity { get; }

        public string Name { get; set; }
        public string Fmt { get; set; }
        public string FmtElt { get; set; }

        public U Cache { get; protected set; }

        HashSet<U> Elts { get; }

        public void Add(U e) => Elts.Add(e);
        public bool Contains(U e) => Elts.Contains(e);
        public int Count => Elts.Count;
        public IEnumerable<U> AllElements() => Elts;

        public bool Equals(IFSet other) => GetHashCode() == other.GetHashCode();
    }

    public interface ISubSet<U> : IFSet<U> where U : struct, IElt
    {
        FSet<U> UpperSet { get; }
    }

    public class Eq<U> : EqualityComparer<U> where U : struct, IElt
    {
        public override bool Equals(U x, U y) => x.Equals(y);
        public override int GetHashCode(U obj) => obj.HashCode;
    }

}
