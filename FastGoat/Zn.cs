using System;
using System.Collections.Generic;
using System.Linq;

using FastGoat.GroupTheory;
using FastGoat.SetTheory;

namespace FastGoat
{
    public struct ZnElt : IElt
    {
        public int HashCode { get; }
        public int[] Table { get; }
        public IFSet FSet { get; }
        public Zn Zn { get; }

        public ZnElt(Zn zn, int hash)
        {
            Zn = zn;
            FSet = zn;
            Table = new int[zn.Dims.Length];
            HashCode = hash;
        }

        public ZnElt(Zn zn, int[] arr, int hash)
        {
            Zn = zn;
            FSet = zn;
            Table = arr.ToArray();
            HashCode = hash;
        }

        public int CompareTo(IElt other) => Helpers.ArrayCompare(Table, other.Table);

        public string Display(string name = "")
        {
            var fmt = string.IsNullOrEmpty(name) ? "{0}" : "{1} = {0}";
            return string.Format(fmt, this, name);
        }

        public bool Equals(IElt other) => HashCode == other.HashCode;
        public override int GetHashCode() => HashCode;
        public override string ToString() => string.Format(FSet.FmtElt, EltStr());

        public string EltStr() => Table.Glue(sep: ",");

        public string Infos() => "";
    }

    public class Zn : Group<ZnElt>
    {
        public int[] Dims { get; }

        public Zn(params int[] dims) : base(dims.Prod())
        {
            Dims = dims;
            Cache = new ZnElt(this, 0);
            Name = string.Join(" x ", Dims.Select(n => $"Z/{n}Z"));
            Fmt = "|{0}| = {1} {2}";
            FmtElt = "({0})";
        }

        public override ZnElt Neutral => new ZnElt(this, 0);
        public override ZnElt Invert(ZnElt a)
        {
            int hash = Helpers.InvertModulo(Dims, a.Table, Cache.Table);
            return new ZnElt(this, Cache.Table, hash);
        }

        public override ZnElt Op(ZnElt a, ZnElt b)
        {
            int hash = Helpers.AddModulo(Dims, a.Table, b.Table, Cache.Table);
            return new ZnElt(this, Cache.Table, hash);
        }

        public ZnElt CreateElement(params int[] vs)
        {
            int hash = Helpers.AddModulo(Dims, vs, Neutral.Table, Cache.Table);
            return new ZnElt(this, Cache.Table, hash);
        }
    }
}
