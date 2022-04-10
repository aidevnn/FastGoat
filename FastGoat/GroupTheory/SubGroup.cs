using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.SetTheory;

namespace FastGoat.GroupTheory
{
    public abstract class SubGroup<U> : SubSet<U>, ISubGroup<U> where U : struct, IElt
    {
        protected SubGroup(Group<U> group, string name, string fmt) : base(group, name, fmt)
        {
            UpperGroup = group;
            ElementOrder = new Dictionary<U, int>();
            SortBy = SortBy.Order;
        }

        Dictionary<U, int> ElementOrder { get; set; }

        public int GetOrder(U e) => ElementOrder[e];

        public Group<U> UpperGroup { get; }

        public abstract U Neutral { get; }
        public abstract U Invert(U a);
        public abstract U Op(U a, U b);
        public XOpLR OpLR { get; set; } = XOpLR.Both;

        public bool IsGroup
        {
            get
            {
                foreach (var e0 in Elts)
                {
                    var inv = Invert(e0);
                    if (!Elts.Contains(inv))
                        return false;
                    foreach (var e1 in Elts)
                    {
                        if (!Elts.Contains(Op(e0, e1)))
                            return false;
                    }
                }

                return true;
            }
        }

        public bool IsCommutative
        {
            get
            {
                foreach (var e0 in Elts)
                {
                    foreach (var e1 in Elts)
                    {
                        var e2 = Op(e0, e1);
                        if (!Elts.Contains(e2))
                            return false;

                        var e3 = Op(e1, e0);
                        if (!Elts.Contains(e3))
                            return false;

                        if (!e2.Equals(e3))
                            return false;
                    }
                }

                return true;
            }
        }

        bool DisplayOrders { get; set; }
        void ComputeOrders()
        {
            List<(U e, int o)> orders = new List<(U, int)>();
            foreach (var e in Elts)
            {
                int ord = 0;
                var acc = Neutral;
                while (ord == 0 || !acc.Equals(Neutral))
                {
                    ++ord;
                    if (OpLR == XOpLR.Left)
                    {
                        acc = Op(e, acc);
                        if (!Elts.Contains(acc))
                            return;
                    }
                    else
                    {
                        acc = Op(acc, e);
                        if (!Elts.Contains(acc))
                            return;
                    }
                }

                orders.Add((e, ord));
            }

            ElementOrder = orders.ToDictionary(a => a.e, b => b.o);
        }

        public void Details(bool displayOrders = true)
        {
            DisplayOrders = displayOrders;
            DisplayElements(displayOrders);
            Table();
            Console.WriteLine();
        }

        public void DisplayElements(bool displayOrders = true)
        {
            DisplayOrders = displayOrders;
            base.DisplayElements();
        }

        public override void DisplayHead()
        {
            ComputeOrders();
            base.DisplayHead();
            Console.WriteLine("IsGroup      :{0,6}", IsGroup);
            Console.WriteLine("IsCommutative:{0,6}", IsCommutative);
            Console.WriteLine();
            SkipFirst = !IsGroup;
        }

        public override int EltCompare(U a, U b)
        {
            if(SortBy == SortBy.Order && DisplayOrders && ElementOrder.ContainsKey(a) && ElementOrder.ContainsKey(b))
            {
                var ordA = ElementOrder[a];
                var ordB = ElementOrder[b];
                if (ordA != ordB) return ordA.CompareTo(ordB);
            }

            return base.EltCompare(a, b);
        }

        protected override string DisplayElement(U e, string name)
        {
            if (!DisplayOrders || ElementOrder.Values.All(a => a == 1))
                return string.Format("{0} = ({1}){2}", name, e.EltStr(), e.Infos());

            return string.Format("{0} = ({1})[{2}]{3}", name, e.EltStr(), ElementOrder[e], e.Infos());
        }

        public void Table()
        {
            if (Elts.Count == 0)
            {
                Console.WriteLine("Empty Set");
                return;
            }

            base.DisplayHead();
            if (!IsGroup)
            {
                Console.WriteLine("Not a Group, need to be developed");
                return;
            }

            if (Elts.Count > 50)
            {
                Console.WriteLine("TOO BIG");
                return;
            }

            ComputeOrders();
            var elts = Elts.ToList();
            elts.Sort(EltCompare);

            var word = GenLetters(Elts.Count, SkipFirst).Select(w => w[0]).ToList();
            Dictionary<char, U> ce = new Dictionary<char, U>();
            Dictionary<U, char> ec = new Dictionary<U, char>(new Eq<U>());

            for (int k = 0; k < elts.Count; ++k)
            {
                var c = word[k];
                var e = elts.ElementAt(k);
                ce[c] = e;
                ec[e] = c;
            }

            string MyFormat(string c, string g, List<char> l) => string.Format("{0,2}|{1}", c, string.Join(g, l));

            var head = MyFormat("*", " ", word);
            var line = MyFormat("--", "", Enumerable.Repeat('-', word.Count * 2).ToList());
            Console.WriteLine(head);
            Console.WriteLine(line);

            foreach (var e0 in elts)
            {
                var v0 = ec[e0].ToString();
                var l0 = elts.Select(e1 => ec[Op(e0, e1)]).ToList();
                Console.WriteLine(MyFormat(v0, " ", l0));
            }

            Console.WriteLine();
        }
    }
}
