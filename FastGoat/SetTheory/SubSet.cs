using System;
using System.Collections.Generic;
using System.Linq;

namespace FastGoat.SetTheory
{
    public enum SortBy { Order, Value }

    public abstract class SubSet<U> : ISubSet<U> where U : struct, IElt
    {
        protected static List<string> GenLetters(int n, bool skipFirst = false)
        {
            int skip = skipFirst ? 1 : 0;
            if (n > 50)
                return Enumerable.Range(skip + 1, n).Select(a => $"E{a,2:000}").ToList();

            return "@abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ".Skip(skip).Take(n).Select(c => $"{c}").ToList();
        }

        protected SubSet(FSet<U> fSet, string name, string fmt)
        {
            UpperSet = fSet;
            Elts = new HashSet<U>(fSet.Capacity, new Eq<U>());
            Name = name;
            Fmt = fmt;
            FmtElt = fSet.FmtElt;
        }

        public FSet<U> UpperSet { get; }
        public SortBy SortBy { get; set; } = SortBy.Value;
        protected HashSet<U> Elts { get; set; }
        public void Add(U e)
        {
            if (!e.FSet.Equals(UpperSet)) return;

            UpperSet.Add(e);
            Elts.Add(e);
        }

        public bool Contains(U e) => Elts.Contains(e);
        public int Count => Elts.Count;

        public string Name { get; set; }
        public string Fmt { get; set; }
        public string FmtElt { get; set; }

        public string Infos { get; set; }
        protected bool SkipFirst { get; set; } = false;

        public IEnumerable<U> AllElements() => Elts;
        public bool Equals(IFSet other) => GetHashCode() == other.GetHashCode();

        public virtual int EltCompare(U a, U b) => a.CompareTo(b);

        public virtual void DisplayHead()
        {
            Console.WriteLine(Fmt, Name, Count, Infos);
        }

        protected virtual string DisplayElement(U e, string name) => string.Format("{0} = {1}", name, e);

        public void DisplayElements()
        {
            if (Elts.Count == 0)
            {
                Console.WriteLine("Empty Set");
                return;
            }

            DisplayHead();
            if (Elts.Count > 300)
            {
                Console.WriteLine("TOO BIG");
                return;
            }

            var elts = Elts.ToList();
            elts.Sort(EltCompare);

            var word = GenLetters(elts.Count, SkipFirst);
            for (int k = 0; k < elts.Count; ++k)
                Console.WriteLine(DisplayElement(elts.ElementAt(k), word[k]));

            Console.WriteLine();
        }

        public bool IsEqual(SubSet<U> set)
        {
            if (!UpperSet.Equals(set.UpperSet)) return false;
            return Elts.SetEquals(set.Elts);
        }
    }

    public class EqSubSet<U> : EqualityComparer<SubSet<U>> where U : struct, IElt
    {
        public override bool Equals(SubSet<U> x, SubSet<U> y) => x.IsEqual(y);
        public override int GetHashCode(SubSet<U> obj) => 1;
    }
}
