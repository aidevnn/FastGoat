using System;
using System.Linq;

namespace FastGoat
{
    public class Tuple2Array
    {
        public int[] Table { get; private set; }
        public Tuple2Array(params int[] table)
        {
            Table = table.ToArray();
        }

        public static Tuple2Array Canonic(int n, int rank)
        {
            var arr = new int[n];
            arr[rank] = 1;
            return new Tuple2Array(arr);
        }

        public static Tuple2Array[] BaseCanonic(int n) => n.Range().Select(rk => Canonic(n, rk)).ToArray();

        public static implicit operator Tuple2Array(int a) => new Tuple2Array(a);
        public static implicit operator Tuple2Array((int a, int b) p) => new Tuple2Array(p.a, p.b);
        public static implicit operator Tuple2Array((int a, int b, int c) p) => new Tuple2Array(p.a, p.b, p.c);
        public static implicit operator Tuple2Array((int a, int b, int c, int d) p) => new Tuple2Array(p.a, p.b, p.c, p.d);
        public static implicit operator Tuple2Array((int a, int b, int c, int d, int e) p) => new Tuple2Array(p.a, p.b, p.c, p.d, p.e);
        public static implicit operator Tuple2Array((int a, int b, int c, int d, int e, int f) p) => new Tuple2Array(p.a, p.b, p.c, p.d, p.e, p.f);

        public override string ToString() => string.Format("({0})", string.Join(" ", Table.Select(e => $"{e,2}")));
    }
}
