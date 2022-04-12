using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;

using FastGoat.SetTheory;
using FastGoat.GroupTheory;

namespace FastGoat
{
    class MainClass
    {
        static void SamplesZn()
        {
            var z = new Zn(2, 3, 3);
            var i0 = z.CreateElement(1, 0, 0);
            var i1 = z.CE(0, 1, 0);
            var i2 = z.CE(0, 0, 1);
            z.Union(i0, i1, i2).Develop("G");

            var j0 = z.Monogenic(i0, "H");
            var j1 = z.Monogenic(z.CE(0, 1, 1), "K");
            j0.DirectProduct(j1).Details();
        }

        static void SamplesSn()
        {
            var s = new Sn(4);
            var a = s.Cycles((1, 3), (2, 4));
            var b = s.Cycles((1, 2), (3, 4));
            var c = s.KCycle(4);

            var g = s.Monogenic(a, "G");
            var h = s.Monogenic(b, "H");
            g.DisplayElements();
            h.DisplayElements();
            var k = g.Union(h).SubGroupOf(s);
            k.Details();
            k.Develop("<GuH>").Details();

            s.Union(a, b, c).Develop("K2").Details();
        }

        static void SamplesZnQuotient()
        {
            var z = new Zn(4, 5);
            var g = z.Union(z.CE(1, 0),z.CE(0, 1)).Develop("G");
            var h = z.Monogenic(z.CE(0, 1), "H");
            var gh = g.Over(h);
            g.SortBy = h.SortBy = gh.SortBy = SortBy.Value;
            g.DisplayElements();
            h.DisplayElements();
            gh.Details();
            gh.DisplayClasses();
        }

        static void SamplesSnQuotient()
        {
            var S4 = new Sn(4);
            var A4 = S4.Union(S4.Cycle(1, 2, 3), S4.Cycle(2, 3, 4)).Develop("A4");
            var Klein = S4.Union(S4.Cycles((1, 2), (3, 4)), S4.Cycles((1, 3), (2, 4))).Develop("Klein");
            var Q = A4.Over(Klein);
            A4.DisplayElements();
            Klein.DisplayElements();
            Q.Details();
            Q.DisplayClasses();
        }

        public static void Main(string[] args)
        {
            //SamplesZn();
            //SamplesSn();
            //SamplesZnQuotient();
            SamplesSnQuotient();
        }
    }
}
