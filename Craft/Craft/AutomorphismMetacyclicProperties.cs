using System.Collections;
using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Perms;

namespace Craft.Craft;

public record AutomorphismMetacyclicProperties(
    (int m, int n, int r) coefs,
    bool mnCoprimes,
    bool mnReductibles,
    bool mnReductionCoprimes,
    (string shortName, string fullName) names,
    (Perm a, Perm b) mtGens,
    ConcreteGroup<Perm> MtInit,
    (Perm a, Perm b) mtRegGens,
    ConcreteGroup<Perm> MtReg,
    ConcreteGroup<Perm> AutMt,
    ConcreteGroup<Perm> InvA,
    ConcreteGroup<Perm> AutCn,
    ConcreteGroup<Perm> InvB,
    ConcreteGroup<Perm> AutCm,
    ConcreteGroup<Perm> ProdInvs,
    ConcreteGroup<Perm> ProdAutCmInvA,
    ConcreteGroup<Perm> ProdAutCnInvB,
    bool ProdInvsNormal,
    bool InvANormal,
    bool InvBNormal,
    bool AutCmNormal,
    bool AutCnNormal,
    bool AutCmInvAIsNormal)
{
    public int m => coefs.m;
    public int n => coefs.n;
    public int r => coefs.r;
    public string shortName => names.shortName;
    public string fullName => names.fullName;
    public string Name => string.Equals(shortName, fullName) ? shortName : $"{shortName} = {fullName}";
    public bool InvAisCyclic => InvA.Order == InvA.ElementsOrdersList().Max();
    public bool InvAisAbelian => InvA.GroupType == GroupType.AbelianGroup;
    public bool InvBisCyclic => InvB.Order == InvB.ElementsOrdersList().Max();
    public bool InvBisAbelian => InvB.GroupType == GroupType.AbelianGroup;
    public int QuoAutMtByProdInvs => AutMt.Order / ProdInvs.Order;

    public string InfosNormals
    {
        get
        {
            var fmt = new[] { "ProdInvsNormal", "InvANormal", "InvBNormal", "AutCmNormal", "AutCnNormal" }
                .Select((s, i) => $"{s}:{{{i},-6}}").Glue(" ");
            return string.Format(fmt, ProdInvsNormal, InvANormal, InvBNormal, AutCmNormal, AutCnNormal);
        }
    }

    public string InfosCoefs
    {
        get
        {
            var fmt = new[] { "mnCoprimes", "mnReductibles", "mnReductionCoprimes" }
                .Select((s, i) => $"{s}:{{{i},-6}}").Glue(" ");
            return string.Format(fmt, mnCoprimes, mnReductibles, mnReductionCoprimes);
        }
    }

    public string NameInfos => $"{Name,-40} {InfosCoefs} {InfosNormals}";

    public override string ToString() => Name;
}

public readonly struct MtElt : IElt<MtElt>
{
    public int I { get; }
    public int J { get; }
    public (int I, int J) IJ => (I, J);

    public MtElt()
    {
        (I, J) = (0, 0);
        Hash = (0, 0).GetHashCode();
    }

    public MtElt(int i, int j)
    {
        (I, J) = (i, j);
        Hash = (i, j).GetHashCode();
    }
    public bool Equals(MtElt other) => IJ.Equals(other.IJ);

    public int CompareTo(MtElt other) => IJ.CompareTo(other.IJ);

    public int Hash { get; }
    public override int GetHashCode() => Hash;

    public override string ToString()
    {
        var ai = I == 1 ? "a" : $"a^{I}";
        var bi = J == 1 ? "b" : $"b^{J}";
        return I == 0 && J == 0 ? "1"
            : J == 0 ? ai
            : I == 0 ? bi
            : $"{ai} * {bi}";
    }
}

public readonly struct MtGroup : IGroup<MtElt>
{
    public int M { get; }
    public int N { get; }
    public int R { get; }
    public (int M, int N, int R) MNR => (M, N, R);

    public MtGroup(int m, int n, int r)
    {
        (M, N, R) = (m, n, r);
        Hash = (m, n, r).GetHashCode();
        Name = FG.MetaCyclicName(m, n, r);
    }

    public IEnumerator<MtElt> GetEnumerator() => GetElements().GetEnumerator();

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Equals(IGroup<MtElt>? other) => other?.Hash == Hash;

    public int Hash { get; }
    public string Name { get; }

    public MtElt this[params ValueType[] us] =>
        new(IntExt.AmodP(Convert.ToInt32(us[0]), M), IntExt.AmodP(Convert.ToInt32(us[1]), N));

    public IEnumerable<MtElt> GetElements() => M.Range().Grid2D(N.Range()).Select(e => new MtElt(e.t1, e.t2));

    public IEnumerable<MtElt> GetGenerators()
    {
        yield return new(1, 0);
        yield return new(0, 1);
    }

    public MtElt Neutral() => new();

    public MtElt Invert(MtElt e)
    {
        ZnInt r = new ZnInt(M, R);
        var i = (-e.I) * r.Pow(-e.J);
        ZnInt j = new(N, -e.J);
        return new(i.K, j.K);
    }

    public MtElt Op(MtElt e1, MtElt e2)
    {
        ZnInt r = new(M, R);
        var i = e1.I + e2.I * r.Pow(e1.J);
        ZnInt j = new(N, e1.J + e2.J);
        return new(i.K, j.K);
    }

    public override int GetHashCode() => Hash;
    public override string ToString() => Name;
    public static ConcreteGroup<MtElt> Create(int m, int n, int r) => Group.Create(new MtGroup(m, n, r));
}