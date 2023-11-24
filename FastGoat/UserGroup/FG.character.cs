using System.Numerics;
using FastGoat.Commons;
using FastGoat.Examples;
using FastGoat.Structures;
using FastGoat.Structures.CartesianProduct;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.VecSpace;
using FastGoat.UserGroup.Characters;
using FastGoat.UserGroup.Integers;
using FastGoat.UserGroup.Polynoms;

namespace FastGoat.UserGroup;

public static partial class FG
{
    public static Complex Substitute(this KPoly<Rational> f, Complex c) => f.Coefs.Select((k, i) => k * Complex.Pow(c, i))
        .Aggregate(Complex.Zero, (sum, ci) => sum + ci);

    public static Character<T> CharacterOne<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        return new Character<T>(gr);
    }

    public static CharacterTable<T> CharacterTable<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        return new CharacterTable<T>(gr);
    }

    public static CharacterTable<T> CharacterTableEmpty<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        return new CharacterTable<T>(gr, true);
    }

    public static List<Character<T>> LinearCharacters<T>(ConcreteGroup<T> gr) where T : struct, IElt<T>
    {
        if (gr.GroupType == GroupType.NonAbelianGroup)
            throw new GroupException(GroupExceptionType.GroupDef);

        var facts = AbelianInvariantsFactors.Reduce(gr).ToArray();
        var gr0 = Abelian(facts);
        var isoMap = Group.AllMorphisms(gr, gr0, Group.MorphismType.Isomorphism).First();
        var o = gr.Count();
        var w = new Cnf(o);
        var wi = facts.Select(i => w.Pow(o / i)).ToArray();
        var clG = Group.ConjugacyClasses(gr);
        var rg = facts.Length.Range();

        var ct = new List<Character<T>>();
        foreach (var h1 in gr0)
        {
            var mapChi = clG.ToDictionary(e => e, _ => (Cnf?)Cnf.CnfZero);
            foreach (var g2 in gr)
            {
                var h2 = isoMap[g2];
                mapChi[g2] = rg.Aggregate(w.One, (prod, i) => prod * wi[i].Pow(h1.Ei[i].K * h2.Ei[i].K));
            }

            ct.Add(new(clG, mapChi));
        }

        return ct;
    }

    public static Cnf InnerProduct<T>(Character<T> chi1, Character<T> chi2) where T : struct, IElt<T>
    {
        if (!chi1.Gr.Equals(chi2.Gr))
            throw new GroupException(GroupExceptionType.GroupDef);

        return chi1.Gr.Aggregate(Cnf.CnfZero, (sum, g) => sum + chi1[g]!.Value * chi2[g]!.Value.Conj) / chi1.Gr.Count();
    }

    public static Character<T> Lift<T>(this Character<Coset<T>> chi, ConjugacyClasses<T> clG) where T : struct, IElt<T>
    {
        var quo = (Quotient<T>)chi.Gr.BaseGroup;
        var map = clG.ToDictionary(e => e, e => chi[quo.GetRepresentative(e)]);
        return new(clG, map);
    }

    public static Character<T> Restriction<T>(this Character<T> chi, ConjugacyClasses<T> clH) where T : struct, IElt<T>
    {
        if (!clH.Gr.SubSetOf(chi.Gr))
            throw new GroupException(GroupExceptionType.GroupDef);

        return new(clH, clH.GetRepresentatives().ToDictionary(e => e, e => chi.Map[chi.Classes.GetRepresentative(e)]));
    }

    public static Character<T> Induction<T>(this Character<T> chi, ConjugacyClasses<T> clG) where T : struct, IElt<T>
    {
        if (!chi.Gr.SubSetOf(clG.Gr))
            throw new GroupException(GroupExceptionType.GroupDef);

        var map = clG.ToDictionary(e => e, _ => (Cnf?)Cnf.CnfZero);
        var oH = chi.Gr.Count();
        foreach (var g in clG)
        {
            var sum = Cnf.CnfZero;
            foreach (var y in clG.Gr)
            {
                var yigy = clG.Gr.Op(clG.Gr.Invert(y), clG.Gr.Op(g, y));
                if (chi.Gr.Contains(yigy))
                    sum += chi[yigy]!.Value;
            }

            map[g] = (sum / oH);
        }

        return new(clG, map);
    }
}