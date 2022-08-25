using FastGoat.Structures.GroupTheory;
using FastGoat.Structures.SetTheory;

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
            var h = g.Monogenic(z.CE(0, 1));

            g.SortBy = h.SortBy = SortBy.Value;
            g.DisplayElements("G");
            h.DisplayElements("H");

            var gh = g.Over(h);
            gh.SortBy = SortBy.Value;
            gh.DisplayElements();
            gh.DisplayGroupTable();
            gh.DisplayClasses();
        }

        static void SamplesZnComputing()
        {
            var z = new Zn(2, 2, 2, 3);
            var z24 = z.GenerateAll();
            z24.DisplayElements("G", "Cartesian product Z/2Z x Z/2Z x Z/2Z x Z/3Z");

            // Greatest order element of the group
            var c6 = z24.Monogenic(z.CE(1, 1, 1, 2));
            c6.DisplayElements("C6");

            // Quotient group 
            var q0 = z24.Over(c6);
            q0.Details();

            // Greatest order element of the quotient group
            var c20 = q0.Monogenic(z.CE(0, 0, 1, 0));
            c20.DisplayElements("C2");

            var q1 = q0.Over(c20);
            q1.Details();

            var c21 = q1.Monogenic(z.CE(0, 1, 0, 0));
            c21.DisplayElements("C2'");

            // Direct product of the elementaries divisors
            Console.WriteLine("###########");
            c6.DirectProduct(c20).DisplayElements("C6.C2");
            Console.WriteLine("###########");
            c6.DirectProduct(c20).DirectProduct(c21).DisplayElements("C6.C2.C2'");
        }

        static void SamplesSn()
        {
            var s = new Sn(4);
            var a = s.C((1, 3), (2, 4));
            var b = s.C((1, 2), (3, 4));
            var c = s.KCycle(4);
            var s0 = s.GroupElement(a, b, c);

            var g = s0.Monogenic(a);
            var h = s0.Monogenic(b);
            g.DisplayElements("G");
            h.DisplayElements("H");
            var k = g.GroupUnion(h);
            k.DisplayElements("GuH");
            k.Generate().Details("<GuH>");

            s.GroupElement(a, b, c).Generate().Details("K2");

            Console.WriteLine($"k:{k.Count}");
            Console.WriteLine($"S0:{s0.Count}");
            Console.WriteLine($"S:{s.Count}");
        }

        static void SamplesSnQuotient()
        {
            var s3 = new Sn(3);
            var S3all = s3.GroupElement(s3.C(1, 2, 3), s3.C(1, 2)).Generate();
            var C3 = S3all.Monogenic(s3.C(1, 2, 3));
            S3all.DisplayElements("S3");
            C3.DisplayElements("C3", "in S3");

            var Q3 = S3all.Over(C3);
            Q3.Details(infos: "in S3");
            Q3.DisplayClasses();

            var S4 = new Sn(4);
            var A4 = S4.GroupElement(S4.C(1, 2, 3), S4.C(2, 3, 4)).Generate();
            var Klein = A4.GroupElement(S4.C((1, 2), (3, 4)), S4.C((1, 3), (2, 4))).Generate();
            A4.DisplayElements("A4", "in S4");
            Klein.DisplayElements("Klein", "in S4");

            var Q4 = A4.Over(Klein);
            Q4.Details(infos: "in S4");
            Q4.DisplayClasses();
        }

        static void SamplesZnInvariants()
        {
            var z14x21 = Zn.CartesianProduct(14, 21);
            GroupExt.InvariantFactors(z14x21);

            var z20x30 = Zn.CartesianProduct(20, 30);
            GroupExt.InvariantFactors(z20x30);

            var z8x18x30 = Zn.CartesianProduct(8, 18, 30); // May the BRUTEFORCE be with you
            GroupExt.InvariantFactors(z8x18x30);
        }

        static void SamplesSnInvariants()
        {
            var s6 = new Sn(6);
            var H0 = s6.GroupElement(s6.KCycle(2), s6.KCycle(3, 4)).Generate("H0");
            H0.DisplayElements();
            GroupExt.InvariantFactors(H0);

            var s9 = new Sn(9);
            var H1 = s9.GroupElement(s9.KCycle(3), s9.KCycle(4, 6)).Generate("H1");
            H1.DisplayElements();
            GroupExt.InvariantFactors(H1);
        }

        public static void Main(string[] args)
        {
            // SamplesZn();
            // SamplesZnQuotient();
            // SamplesZnComputing();
            // SamplesZnInvariants();

            // SamplesSn();
            // SamplesSnQuotient();
            SamplesSnInvariants();
        }
    }
}
