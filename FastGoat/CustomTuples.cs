using System;
using System.Linq;

namespace FastGoat
{

    public class SingleTuple
    {
        public int[] Table { get; private set; }
        public SingleTuple(params int[] table)
        {
            Table = table.ToArray();
        }

        public static SingleTuple Canonic(int n, int rank)
        {
            var arr = new int[n];
            arr[rank] = 1;
            return new SingleTuple(arr);
        }

        public static SingleTuple[] BaseCanonic(int n) => n.Range().Select(rk => Canonic(n, rk)).ToArray();

        public static implicit operator SingleTuple(int a) => new SingleTuple(a);
        public static implicit operator SingleTuple((int a, int b) p) => new SingleTuple(p.a, p.b);
        public static implicit operator SingleTuple((int a, int b, int c) p) => new SingleTuple(p.a, p.b, p.c);
        public static implicit operator SingleTuple((int a, int b, int c, int d) p) => new SingleTuple(p.a, p.b, p.c, p.d);
        public static implicit operator SingleTuple((int a, int b, int c, int d, int e) p) => new SingleTuple(p.a, p.b, p.c, p.d, p.e);
        public static implicit operator SingleTuple((int a, int b, int c, int d, int e, int f) p) => new SingleTuple(p.a, p.b, p.c, p.d, p.e, p.f);

        public override string ToString() => string.Format("({0})", string.Join(" ", Table.Select(e => $"{e,2}")));
    }

    public class ManyTuples
    {
        public SingleTuple[] Tuples { get; private set; }
        public ManyTuples(params SingleTuple[] tuples)
        {
            Tuples = tuples.ToArray();
        }

        public static implicit operator ManyTuples(int a) => new ManyTuples(a);
        public static implicit operator ManyTuples((int a, int b) p) => new ManyTuples(p);
        public static implicit operator ManyTuples((int a, int b, int c) p) => new ManyTuples(p);
        public static implicit operator ManyTuples((int a, int b, int c, int d) p) => new ManyTuples(p);
        public static implicit operator ManyTuples((int a, int b, int c, int d, int e) p) => new ManyTuples(p);
        public static implicit operator ManyTuples((int a, int b, int c, int d, int e, int f) p) => new ManyTuples(p);

        public static implicit operator ManyTuples(SingleTuple mc) => new ManyTuples(mc);
        public static implicit operator ManyTuples((SingleTuple a, SingleTuple b) mc) => new ManyTuples(mc.a, mc.b);
        public static implicit operator ManyTuples((SingleTuple a, SingleTuple b, SingleTuple c) mc) => new ManyTuples(mc.a, mc.b, mc.c);
        public static implicit operator ManyTuples((SingleTuple a, SingleTuple b, SingleTuple c, SingleTuple d) mc) => new ManyTuples(mc.a, mc.b, mc.c, mc.d);
        public static implicit operator ManyTuples((SingleTuple a, SingleTuple b, SingleTuple c, SingleTuple d, SingleTuple e) mc) => new ManyTuples(mc.a, mc.b, mc.c, mc.d, mc.e);
        public static implicit operator ManyTuples((SingleTuple a, SingleTuple b, SingleTuple c, SingleTuple d, SingleTuple e, SingleTuple f) mc) => new ManyTuples(mc.a, mc.b, mc.c, mc.d, mc.e, mc.f);

        public override string ToString() => string.Format("[{0}]", string.Join(" ", Tuples.Select(e => $"{e}")));
    }
}
