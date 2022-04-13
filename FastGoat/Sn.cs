using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.GroupTheory;
using FastGoat.SetTheory;

namespace FastGoat
{
    public struct Permutation : IElt
    {
        public int HashCode { get; }
        public int[] Table { get; }
        public IFSet FSet { get; }
        public Sn Sn { get; }

        public Permutation(Sn sn)
        {
            Sn = sn;
            FSet = sn;
            Table = sn.N.Range();
            HashCode = sn.NHash;
            Invariants = Table.ToArray();
            Sgn = 1;
            Cycles = new int[0][];
        }

        public Permutation(Sn sn, int[] arr, int hash)
        {
            Sn = sn;
            FSet = sn;
            Table = arr.ToArray();
            HashCode = hash;
            var orbits = Table.Orbits();
            Invariants = orbits.Where(a => a.Count == 1).SelectMany(a => a).ToArray();
            Cycles = orbits.Where(a => a.Count > 1).Select(a => a.ToArray()).ToArray();
            Sgn = (int)Math.Pow(-1, Sn.N - orbits.Count);
        }

        public int Sgn { get; }
        public int[] Invariants { get; }
        public int[][] Cycles { get; }

        public int CompareTo(IElt other) => Helpers.ArrayCompare(Table, other.Table);

        public string Display(string name = "")
        {
            var fmt = string.IsNullOrEmpty(name) ? "{0}" : "{1} = {0}";
            return string.Format(fmt, this, name);
        }

        public bool Equals(IElt other) => HashCode == other.HashCode;
        public override int GetHashCode() => HashCode;
        public override string ToString() => string.Format(FSet.FmtElt, EltStr());

        public string EltStr() => Table.Select(a => 1 + a).Glue(sep: " ");

        public string Infos()
        {
            var inv = Invariants.Length == 0 ? "" : string.Format(" Invariants : [{0}]", Invariants.Select(b => b + 1).Glue(fmt: "{0}"));
            var strCycles = Cycles.Select(a => string.Format("({0})", a.Select(b => b + 1).Glue(fmt: "{0}"))).Glue();
            var cycles = Cycles.Length == 0 ? "" : string.Format(" Cycles : {0}", strCycles);
            var sgn = Sgn == 1 ? "(+)" : "(-)";
            return string.Format("{0}{1}{2}", sgn, inv, cycles);
        }

    }
    public class Sn : Group<Permutation>
    {
        public int N { get; }
        public int NHash { get; }

        public Sn(int n) : base(n.FactorielN())
        {
            N = n;
            NHash = Helpers.GenHash(n, n.Range());

            Cache = new Permutation(this);
            Name = $"S{N}";
            Fmt = "|{0}| = {1} {2}";
            FmtElt = "({0})";
        }

        public override Permutation Neutral => new Permutation(this);

        public override Permutation Invert(Permutation a)
        {
            var hash = Helpers.InvertPermutation(a.Table, Cache.Table);
            return new Permutation(this, Cache.Table, hash);
        }

        public override Permutation Op(Permutation a, Permutation b)
        {
            var hash = Helpers.ComposePermutation(a.Table, b.Table, Cache.Table);
            return new Permutation(this, Cache.Table, hash);
        }

        public Permutation CreateElement(params int[] vs)
        {
            vs.Add(-1);
            if (!Helpers.CheckTable(N, vs))
                return Neutral;

            var hash = Helpers.GenHash(N, vs);
            return new Permutation(this, vs, hash);
        }

        public Permutation Cycle(params int[] vs)
        {
            vs.Add(-1);
            if (!Helpers.CheckCycle(N, vs))
                return Neutral;

            Neutral.Table.CopyTo(Cache.Table, 0);
            Helpers.ApplyCycle(Cache.Table, vs);
            var hash = Helpers.GenHash(N, Cache.Table);
            return new Permutation(this, Cache.Table, hash);
        }

        public Permutation KCycle(int count) => Cycle(count.Range(1));
        public Permutation KCycle(int start, int count) => Cycle(count.Range(start));

        public Permutation Cycle(Tuple2Array cycle) => Cycle(cycle.Table);
        public Permutation ComposesCycles(params Tuple2Array[] cycles)
        {
            Neutral.Table.CopyTo(Cache.Table, 0);
            foreach(var e in cycles)
            {
                e.Table.Add(-1);
                if (!Helpers.CheckCycle(N, e.Table)) return Neutral;

                Helpers.ApplyCycle(Cache.Table, e.Table);
            }

            var hash = Helpers.GenHash(N, Cache.Table);
            return new Permutation(this, Cache.Table, hash);
        }

        public Permutation C(params int[] vs) => Cycle(vs);
        public Permutation C(params Tuple2Array[] cycles) => ComposesCycles(cycles);

        public SubGroup<Permutation> GenerateAll()
        {
            return this.GeneratedSubGroup(Helpers.AllPerms(N, true).Select(CreateElement).ToArray()).GName($"S{N}");
        }
    }
}
