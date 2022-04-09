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
        public static void Main(string[] args)
        {
            var s = new Sn(4);
            var a = s.Cycle((1, 3), (2, 4));
            var b = s.Cycle((1, 2), (3, 4));
            var c = s.KCycle(4);

            var g = s.Monogenic(a, "G");
            var h = s.Monogenic(b, "H");
            g.Details();
            h.Details();
            var k = g.Union(h).SubGroupOf(s);
            k.Details();
            k.DirectProduct(k, "K").Details();
            k.Develop("K1").Details();

            s.Union(a, b, c).Develop("K2").Details();
        }
    }
}
