using FastGoat.Commons;

namespace FastGoat.Structures.GenericGroup;

public class IsomorphEquality<T> : EqualityComparer<ConcreteGroup<T>> where T : struct, IElt<T>
{
    public override bool Equals(ConcreteGroup<T>? x, ConcreteGroup<T>? y)
    {
        return x is not null && y is not null && x.IsIsomorphicTo(y);
    }

    public override int GetHashCode(ConcreteGroup<T> obj) => (obj.Count(), obj.GroupType).GetHashCode();
}

public class IsomorphSubGroupsInfosEquality<T> : EqualityComparer<AllSubgroups<T>> where T : struct, IElt<T>
{
    public override bool Equals(AllSubgroups<T> x, AllSubgroups<T> y)
    {
        return x.Parent.IsIsomorphicTo(y.Parent);
    }

    public override int GetHashCode(AllSubgroups<T> obj) => obj.Infos.GetHashCode();
}

public class OpByAutEquality<T1, T2>(ConcreteGroup<T2> G, ConcreteGroup<Automorphism<T2>> AutG, ConcreteGroup<Automorphism<T1>> AutN)
    : EqualityComparer<Homomorphism<T2, Automorphism<T1>>>
    where T1 : struct, IElt<T1>
    where T2 : struct, IElt<T2>
{
    public override bool Equals(Homomorphism<T2, Automorphism<T1>> x, Homomorphism<T2, Automorphism<T1>> y)
    {
        var N = ((AutomorphismGroup<T1>)AutN.BaseGroup).G;
        var GN = G.Grid2D(N).Select(e => (g: e.t1, n: e.t2)).ToArray();

        return AutG.Any(f => GN.All(e => x[e.g][e.n].Equals(y[f[e.g]][e.n])))
               || AutN.Any(f =>
               {
                   var fi = f.Invert();
                   return GN.All(e => f[x[e.g][fi[e.n]]].Equals(y[e.g][e.n]));
               });
    }

    public override int GetHashCode(Homomorphism<T2, Automorphism<T1>> obj) => (typeof(T1), typeof(T2), obj.Count).GetHashCode();
}
