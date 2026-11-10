using FastGoat.Commons;
using FastGoat.Structures;
using FastGoat.Structures.GenericGroup;
using FastGoat.Structures.Naming;
using FastGoat.Structures.Subgroups;

namespace Craft.Craft;

public struct PolyCyclicSequence<T>(List<(ConcreteGroup<Coset<T>> K, T e, bool Split)> gens) where T : struct, IElt<T>
{
    public List<ConcreteGroup<T>> Normals => gens.Select(e => e.K.Neutral().H).ToList();
    public List<(ConcreteGroup<Coset<T>> K, T e, bool Split)> Generators => gens;
    public PolyCyclicSequence<T> Add(ConcreteGroup<Coset<T>> K, T e, bool Split) => new(Generators.Prepend((K, e, Split)).ToList());

    public void FactorElements()
    {
        var G = Normals.Last();
        var all = Generators.SkipLast(1).Select(e => e.K.Order.SeqLazy().ToDictionary(i => i, _ => e.e).ToArray()).MultiLoop()
            .Select(l => l.ToArray()).ToArray();
        all.Println(l =>
            $"{l.Select(e => $"({e.Value} ^ {e.Key})").Glue(" * "),-50} => {G.OpSeq(l.Select(e => G.Times(e.Value, e.Key)))}");
    }

    public string Relators()
    {
        var G = Normals.Last();
        var gens0 = Generators.SkipLast(1).Reverse().ToArray();
        var allNormalForms = gens0.Select(e => e.K.Order.SeqLazy().ToDictionary(i => i, _ => e.e).ToArray()).MultiLoop()
            .ToDictionary(l => G.OpSeq(l.Select(e => G.Times(e.Value, e.Key))), l => l.ToArray());
        var Xis = gens0.Length.SeqLazy().ToDictionary(i => i, i => $"{(char)('a' + i)}");

        string Format(int k, KeyValuePair<int, T> kv)
        {
            var xi = Xis[k];
            var j = kv.Key;
            return j == 0 ? "" : j == 1 ? xi : $"{xi}{j}";
        }

        var relators = new Dictionary<string, string>();
        foreach (var i in gens0.Length.SeqLazy())
        {
            var xi = gens0[i].e;
            var ord = G.ElementsOrders[xi];
            
            var Rii = ord.SeqLazy(1).Select(j => (j, xi_j: G.Times(xi, j))).Select(e => (e.j, e.xi_j, seq: allNormalForms[e.xi_j]))
                .First(e => e.seq.Take(i + 1).All(f => f.Key == 0));
            var strRii = $"{Rii.seq.Index().Select(e => Format(e.Index, e.Item)).Glue()}";
            var rii = string.IsNullOrEmpty(strRii) ? "" : $" = {strRii}";
            // Console.WriteLine($"R{i}{i}: {Xis[i]}{Rii.j}{rii}");
            relators[$"R{i}{i}"] = $"{Xis[i]}{Rii.j}{rii}";
        }

        foreach (var i in gens0.Length.SeqLazy().Skip(1))
        {
            var xi = gens0[i].e;
            foreach (var Rij in gens0.Length.SeqLazy().Where(j => j < i)
                         .Select(j => (j, xj1_xi_xj: G.OpSeq(G.Invert(gens0[j].e), xi, gens0[j].e)))
                         .Select(e => (e.j, e.xj1_xi_xj, seq: allNormalForms[e.xj1_xi_xj]))
                         .Where(e => e.seq.Take(e.j + 1).All(f => f.Key == 0)))
            {
                var strRij =
                    $"{Rij.seq.Index().Select(e => Format(e.Index, e.Item)).Glue()}";
                var rij = string.IsNullOrEmpty(strRij) ? "" : $" = {strRij}";
                // Console.WriteLine($"R{i}{Rij.j}: {Xis[Rij.j]}-1{Xis[i]}{Xis[Rij.j]}{rij}");
                relators[$"R{i}{Rij.j}"] = $"{Xis[Rij.j]}-1{Xis[i]}{Xis[Rij.j]}{rij}";
            }
        }

        foreach (var i in gens0.Length.SeqLazy().Skip(1))
        {
            var xi = gens0[i].e;
            foreach (var Rji in gens0.Length.SeqLazy().Where(j => j < i)
                         .Select(j => (j, xj_xi_xj1: G.OpSeq(gens0[j].e, xi, G.Invert(gens0[j].e))))
                         .Select(e => (e.j, xj1_xi_xj: e.xj_xi_xj1, seq: allNormalForms[e.xj_xi_xj1]))
                         .Where(e => e.seq.Take(e.j + 1).All(f => f.Key == 0)))
            {
                var strRji = $"{Rji.seq.Index().Select(e => Format(e.Index, e.Item)).Glue()}";
                var rji = string.IsNullOrEmpty(strRji) ? "" : $" = {strRji}";
                // Console.WriteLine($"R{Rji.j}{i}: {Xis[Rji.j]}{Xis[i]}{Xis[Rji.j]}-1{rji}");
                // relators[$"R{Rji.j}{i}"] = $"{Xis[Rji.j]}{Xis[i]}{Xis[Rji.j]}-1{rji}";
            }
        }

        Console.WriteLine();
        var strRelators = relators.Values.Distinct().Glue(", ");
        Console.WriteLine("PCGroup Relators: {0}", relators.Values.Distinct().Glue(", "));
        Console.WriteLine();
        return strRelators;
    }

    public static IEnumerable<XSet<T>> KerHom(ConcreteGroup<T> G, ConcreteGroup<T> H)
    {
        return GroupCraft.AllMorphismsWithPruning(G, H)
            .Where(hom => !hom.IsNull && H.Order == hom.Image().Count())
            .Select(hom => hom.Kernel().ToXSet());
    }

    public static IEnumerable<(ConcreteGroup<Coset<T>> K, T e, bool Split)> DecompStep(ConcreteGroup<T> G, HashSet<T> K)
    {
        if (G.GroupType == GroupType.AbelianGroup)
        {
            var e = G.MaxBy(e => G.ElementsOrders[e]);
            var E = Group.Generate("E", G, e);
            var N = Group.FactorGroup(G, E);
            var F = G.Over2(N);
            return [(F, e, true)];
        }
        else
            return G.ElementsOrders.Where(e => IntExt.PrimesDec(e.Value).Count == 1 && !K.Contains(e.Key))
                .OrderBy(e => e.Key).ThenBy(e => e.Value).Select(e => e.Key)
                .DistinctBy(e => Group.Orbits(G.GetGenerators().ToHashSet(), Group.ByConjugate(G), e), new SetEquality<T>())
                .DistinctBy(e => Group.GenerateElements(G, e), new SetEquality<T>())
                .Select(e => (K: Group.Generate("K", G, e), e))
                .SelectMany(e => KerHom(G, e.K).Select(N => (N, e.K, e.e)))
                .Select(e => (e.N, e.K, e.e, inter: e.N.Intersect(e.K).Count()))
                .OrderBy(e => (e.inter, e.K.Order))
                .DistinctBy(e => e.N)
                .Select(e =>
                {
                    var N = Group.Generate("N", G, e.N.ToArray());
                    var F = G.Over2(N);
                    var x = F.GetGenerators().First().X;
                    return (F, x, e.inter == 1);
                });
    }

    public static IEnumerable<PolyCyclicSequence<T>> Decomp(ConcreteGroup<T> G)
    {
        var all = new List<(ConcreteGroup<T> G, PolyCyclicSequence<T> Partial)>(10 * G.Order) 
            { (G, new([(G.Over2(G), G.Neutral(), true)])) };
        while (all.Count != 0)
        {
            var tmp = all.OrderBy(pol => pol.Partial.Generators.Count(e => !e.Split))
                .ThenBy(pol => pol.Partial.Generators.Count)
                .ThenByDescending(pol => pol.G.Order)
                .ToArray();
            all.Clear();

            foreach (var (Gi, partialDec) in tmp)
            {
                var gensKi = partialDec.Generators.Select(e => e.e).Append(G.Neutral()).ToArray();
                var Ki = GroupCraft.SubGroupsConjugates(G, new GroupSubset<T>([..gensKi], Group.GenerateElements(G, gensKi)))
                    .SelectMany(k => k).ToHashSet();
                foreach (var (K,e,split) in DecompStep(Gi, Ki))
                {
                    var N = K.Neutral().H;
                    var dec = partialDec.Add(K, e, split);
                    if (N.Order == 1)
                        yield return dec;
                    else
                        all.Add((N, dec));
                }
            }
        }
    }
}