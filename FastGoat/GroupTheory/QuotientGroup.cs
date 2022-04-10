using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public class QuotientGroup<U> : SubGroup<U> where U : struct, IElt
    {
        public QuotientGroup(SubGroup<U> g, SubGroup<U> h, XOpLR opLR, string name) : base(g.UpperGroup, name, g.UpperGroup.Fmt)
        {
            G = g;
            H = h;
            OpLR = opLR;
            representatives = new Dictionary<U, U>(g.Count, new Eq<U>());
            classOf = new Dictionary<U, List<U>>(g.Count / h.Count, new Eq<U>());

            Init();
        }

        public QuotientGroup(SubGroup<U> g, SubGroup<U> h, XOpLR opLR) : this(g, h, opLR, $"{g.Name}/{h.Name}") { }
        public QuotientGroup(SubGroup<U> g, SubGroup<U> h) : this(g, h, XOpLR.Both, $"{g.Name}/{h.Name}") { }

        readonly Dictionary<U, U> representatives;
        readonly Dictionary<U, List<U>> classOf;
        SubGroup<U> G { get; }
        SubGroup<U> H { get; }

        void Init()
        {
            Fmt = string.Format("|{{0}}| = {{1}} with |{0}| = {2} and |{1}| = {3}, Op{4}", G.Name, H.Name, G.Count, H.Count, OpLR);

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
            foreach (var e0 in G.AllElements())
            {
                foreach (var e1 in G.AllElements())
                {
                    var e2 = UpperGroup.Op(e0, UpperGroup.Invert(e1));
                    if (H.Contains(e2))
                        equiv.Add((e0, e1));
                }
            }

            foreach (var (x, y) in equiv)
            {
                foreach(var a in G.AllElements())
                {
                    if (OpLR == XOpLR.Left)
                    {
                        var ax = UpperGroup.Op(a, x);
                        var ay = UpperGroup.Op(a, y);
                        var e2 = UpperGroup.Op(ax, UpperGroup.Invert(ay));
                        if (!H.Contains(e2))
                            return false;
                    }
                    else
                    {
                        var xa = UpperGroup.Op(x, a);
                        var ya = UpperGroup.Op(y, a);
                        var e2 = UpperGroup.Op(xa, UpperGroup.Invert(ya));
                        if (!H.Contains(e2))
                            return false;

                    }
                }
            }

            return true;
        }

        void Classes()
        {
            var comp = OpCompatibility();
            if (!comp)
                return;

            HashSet<SubSet<U>> hs0 = new HashSet<SubSet<U>>(new EqSubSet<U>());
            HashSet<SubSet<U>> hs1 = new HashSet<SubSet<U>>(new EqSubSet<U>());
            var listG = G.AllElements().ToList();
            listG.Sort(G.EltCompare);
            foreach(var e0 in listG)
            {
                if (OpLR == XOpLR.Left || OpLR == XOpLR.Both)
                {
                    var equiv = new GroupOp<U>(e0, H);
                    hs0.Add(equiv);
                }

                if (OpLR == XOpLR.Right || OpLR == XOpLR.Both)
                {
                    var equiv = new GroupOp<U>(H, e0);
                    hs1.Add(equiv);
                }
            }

            if (OpLR == XOpLR.Right)
                hs0 = new HashSet<SubSet<U>>(hs1, new EqSubSet<U>());

            if (OpLR == XOpLR.Both)
            {
                if (hs0.Count != hs1.Count || hs0.Any(lh => !hs1.Contains(lh)))
                    return;
            }

            foreach (var lh in hs0)
            {
                var lu = new List<U>(lh.AllElements());
                lu.Sort(G.EltCompare);
                var r = lu.First();
                var subG = new Union<U>(G.UpperSet, lu.ToArray());
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
