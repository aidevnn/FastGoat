using FastGoat.Structures.GroupTheory;

namespace FastGoat
{
    public class Program
    {
        static void SamplesZn()
        {
            var z = new Zn(2, 3, 3);
            var i0 = z.CreateElement(1, 0, 0);
            var i1 = z.CE(0, 1, 0);
            var i2 = z.CE(0, 0, 1);
            z.GroupElement(i0, i1, i2).Generate();

            var h = z.Monogenic(i0);
            var k = z.Monogenic(z.CE(0, 1, 1));
            var hk = h.DirectProduct(k);

            h.DisplayElements("H");
            k.DisplayElements("K");
            hk.DisplayElements("H.K");
            hk.DisplayGroupTable();
        }

        static void SamplesZnQuotient()
        {
            var z = new Zn(4, 5);
            var g = z.GroupElement(z.CE(1, 0), z.CE(0, 1)).Generate();
            var h = z.Monogenic(z.CE(0, 1));

            g.SortBy = h.SortBy = SortBy.Value;
            g.DisplayElements("G");
            h.DisplayElements("H");

            var gh = g.Over(h);
            gh.SortBy = SortBy.Value;
            gh.DisplayElements();
            gh.DisplayGroupTable();
            gh.DisplayClasses();
        }

        static void SamplesZnInvFactors()
        {
            var z = new Zn(2, 2, 2, 3);
            var z24 = z.GenerateAll();
            z24.DisplayElements("G", "Cartesian product Z/2Z x Z/2Z x Z/2Z x Z/3Z");

            // Greatest order element of the group
            var c6 = z.Monogenic(z.CE(1, 1, 1, 2));
            c6.DisplayElements("C6");

            // Quotient group 
            var k = z24.Over(c6);
            k.Details();

            // Greatest order element of the quotient group
            var c20 = z.Monogenic(z.CE(0, 0, 1, 0));
            c20.DisplayElements("C2");

            k.Over(c20).Details();

            var c21 = z.Monogenic(z.CE(0, 1, 0, 0));
            c21.DisplayElements("C2'");

            // Direct product of the invariants factors
            c20.DirectProduct(c21).DirectProduct(c6).DisplayElements("C2.C2'.C6");
        }

        static void SamplesSn()
        {
            var s = new Sn(4);
            var a = s.C((1, 3), (2, 4));
            var b = s.C((1, 2), (3, 4));
            var c = s.KCycle(4);

            var g = s.Monogenic(a);
            var h = s.Monogenic(b);
            g.DisplayElements("G");
            h.DisplayElements("H");
            var k = s.GroupUnion(g, h);
            k.DisplayElements();
            k.Generate().Details("<GuH>");

            s.GroupElement(a, b, c).Generate().Details("K2");
        }

        static void SamplesSnQuotient()
        {
            var S4 = new Sn(4);
            var A4 = S4.GroupElement(S4.C(1, 2, 3), S4.C(2, 3, 4)).Generate();
            var Klein = S4.GroupElement(S4.C((1, 2), (3, 4)), S4.C((1, 3), (2, 4))).Generate();
            A4.DisplayElements("A4", "in S4");
            Klein.DisplayElements("Klein", "in S4");

            var Q = A4.Over(Klein);
            Q.Details(infos: "in S4");
            Q.DisplayClasses();
        }

        static void Samples()
        {
            // SamplesZn();
            // SamplesZnQuotient();
            // SamplesZnInvFactors();

            // SamplesSn();
            SamplesSnQuotient();
        }
        public static void Main(string[] args)
        {
            SamplesZn();
            SamplesZnQuotient();
            SamplesZnInvFactors();
            SamplesSn();
            SamplesSnQuotient();
        }
    }
}
