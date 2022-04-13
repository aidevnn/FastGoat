using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public class QuotientGroup<U> : SubGroup<U> where U : struct, IElt
    {
        public QuotientGroup(SubGroup<U> g, SubGroup<U> h) : base(g.UpperGroup)
        {
            G = g;
            H = h;
            representatives = new Dictionary<U, U>(g.Count, new Eq<U>());
            classOf = new Dictionary<U, List<U>>(g.Count / h.Count, new Eq<U>());

            Name = $"{g.Name}/{h.Name}";
            Init();
        }

        readonly Dictionary<U, U> representatives;
        readonly Dictionary<U, List<U>> classOf;
        SubGroup<U> G { get; }
        SubGroup<U> H { get; }

        void Init()
        {
            Fmt = string.Format("|{{0}}| = {{1}} with |{0}| = {2} and |{1}| = {3}", G.Name, H.Name, G.Count, H.Count);

            if (!G.UpperGroup.Equals(H.UpperGroup))
                return;

            if (!G.IsGroup || !H.IsGroup)
                return;

            var h = H.AllElements().ToHashSet(new Eq<U>());
            if (!h.IsSubsetOf(G.AllElements()))
                return;

            HashSet<SubSet<U>> xH = new HashSet<SubSet<U>>(new EqSubSet<U>());
            HashSet<SubSet<U>> Hx = new HashSet<SubSet<U>>(new EqSubSet<U>());
            var listG = G.AllElements().ToList();
            listG.Sort(G.EltCompare);
            foreach (var e0 in listG)
            {
                xH.Add(new GroupOp<U>(e0, H));
                Hx.Add(new GroupOp<U>(H, e0));
            }

            if (!xH.SetEquals(Hx))
                return;

            foreach (var lh in xH)
            {
                var lu = new List<U>(lh.AllElements());
                lu.Sort(G.EltCompare);
                var r = lu.First();
                classOf[r] = lu;
                Add(r);
                lu.ForEach(e => representatives[e] = r);
            }
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a)
        {
            var ia = G.Invert(a);
            return representatives[ia];
        }

        public override U Op(U a, U b)
        {
            var c = G.Op(a, b);
            return representatives[c];
        }

        public U ClassOf(U e) => representatives[e];

        public void DisplayClasses()
        {
            foreach(var kp in classOf)
            {
                Console.WriteLine("Class of : {0}", kp.Key);
                foreach (var e in kp.Value)
                    Console.WriteLine("\t{0}", e);
            }

            Console.WriteLine();
        }
    }
}
