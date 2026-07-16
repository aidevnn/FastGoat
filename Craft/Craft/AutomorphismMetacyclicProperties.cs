using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.UserGroup;
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