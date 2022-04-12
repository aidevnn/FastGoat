using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public class QuotientGroup<U> : SubGroup<U> where U : struct, IElt
    {
        public QuotientGroup(SubGroup<U> g, SubGroup<U> h, string name) : base(g.UpperGroup, name, g.UpperGroup.Fmt)
        {
            G = g;
            H = h;
            representatives = new Dictionary<U, U>(g.Count, new Eq<U>());
            classOf = new Dictionary<U, List<U>>(g.Count / h.Count, new Eq<U>());

            Init();
        }

        public QuotientGroup(SubGroup<U> g, SubGroup<U> h) : this(g, h, $"{g.Name}/{h.Name}") { }

        readonly Dictionary<U, U> representatives;
        readonly Dictionary<U, List<U>> classOf;
        SubGroup<U> G { get; }
        SubGroup<U> H { get; }

        void Init()
        {
            Fmt = string.Format("|{{0}}| = {{1}} with |{0}| = {2} and |{1}| = {3}", G.Name, H.Name, G.Count, H.Count);

            var gIsGr = G.IsGroup;
            var hIsGr = H.IsGroup;

            if (!G.UpperGroup.Equals(H.UpperGroup))
                return;

            if (!gIsGr || !hIsGr)
                return;

            foreach (var e in H.AllElements())
            {
                if (!G.Contains(e))
                {
                    Console.WriteLine(e);
                    return;
                }
            }

            Classes();
        }

        bool OpCompatibility()
        {
            List<(U, U)> equiv = new List<(U, U)>();
            var elts = G.AllElements().ToList();
            int n = elts.Count;
            foreach (var e0 in G.AllElements())
            {
                foreach (var e1 in G.AllElements())
                {
                    if (e0.Equals(e1)) continue;

                    var e20 = UpperGroup.Op(UpperGroup.Invert(e0), e1);
                    var e21 = UpperGroup.Op(e0, UpperGroup.Invert(e1));
                    if (H.Contains(e20) && H.Contains(e21))
                        equiv.Add((e0, e1));
                }
            }

            foreach (var (x, y) in equiv)
            {
                foreach(var a in G.AllElements())
                {
                    var ax = UpperGroup.Op(a, x);
                    var ay = UpperGroup.Op(a, y);
                    var el = UpperGroup.Op(UpperGroup.Invert(ax), ay);
                    if (!H.Contains(el))
                        return false;

                    var xa = UpperGroup.Op(x, a);
                    var ya = UpperGroup.Op(y, a);
                    var er = UpperGroup.Op(xa, UpperGroup.Invert(ya));
                    if (!H.Contains(er))
                        return false;
                }
            }

            return true;
        }

        void Classes()
        {
            var comp = OpCompatibility();
            if (!comp)
                return;

            HashSet<SubSet<U>> xH = new HashSet<SubSet<U>>(new EqSubSet<U>());
            HashSet<SubSet<U>> Hx = new HashSet<SubSet<U>>(new EqSubSet<U>());
            var listG = G.AllElements().ToList();
            listG.Sort(G.EltCompare);
            foreach(var e0 in listG)
            {
                var l = new GroupOp<U>(e0, H);
                xH.Add(l);

                var r = new GroupOp<U>(H, e0);
                Hx.Add(r);
            }

            if (xH.Count != Hx.Count || xH.Any(lh => !Hx.Contains(lh)))
                return;

            foreach (var lh in xH)
            {
                var lu = new List<U>(lh.AllElements());
                lu.Sort(G.EltCompare);
                var r = lu.First();
                classOf[r] = lu;
                Add(r);
                foreach (var e in lu)
                    representatives[e] = r;
            }
        }

        public override U Neutral => UpperGroup.Neutral;
        public override U Invert(U a)
        {
            var ra = representatives[a];
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
