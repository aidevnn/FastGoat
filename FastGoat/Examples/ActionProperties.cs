using FastGoat.Gp;
using FastGoat.UserGroup;

namespace FastGoat.Examples;

public static class ActionProperties
{
    public static void CheckPropertiesForSemiDirectProduct<T1, T2>(ConcreteGroup<T1> grN, ConcreteGroup<T2> grG,
        params T1[] rmPoints)
        where T1 : struct, IElt<T1> where T2 : struct, IElt<T2>
    {
        var autN = Group.AutomorphismGroup(grN);
        var opsGautN = Group.AllHomomorphisms(grG, autN);
        List<Dictionary<T2, Automorphism<T1>>> opsTransitives = new();
        List<Dictionary<T2, Automorphism<T1>>> opsFaithful = new();

        // always remove neutral x = id
        HashSet<T1> fixedPointsSet = rmPoints.Append(grN.Neutral()).ToHashSet();

        foreach (var op in opsGautN)
        {
            var kerOp = op.Where(kp => kp.Value.Equals(autN.Neutral())).Select(kp => kp.Key).ToArray();
            if (kerOp.Length == 1)
                opsFaithful.Add(op);

            var orbs = Group.AllOrbits(grG, grN.Except(fixedPointsSet).ToArray(),
                Group.ByAutomorphism(op));
            if (orbs.Count == 1)
                opsTransitives.Add(op);
        }

        if (opsTransitives.Count == 0)
        {
            foreach (var op in opsGautN)
            {
                Console.WriteLine(
                    "Fixing the transitivity by removing points and setX cardinal divides action group order.");
                Console.WriteLine("Group {0} Acting on Group {1}", grG.Name, grN.Name);
                Console.WriteLine(op.Glue());
                Group.DisplayOrbx(grG, grN.Except(fixedPointsSet).ToArray(), Group.ByAutomorphism(op));
            }
        }

        foreach (var op in opsTransitives)
        {
            Console.WriteLine(
                $"############################# Transitives Group Action Over Set {{{grN.Except(fixedPointsSet).Glue(", ")}}}");
            var sdp = Group.SemiDirectProd(grN, op, grG);
            DisplayGroup.HeadSdp(sdp);
            Console.WriteLine(sdp.ElementsOrdersList().Glue(", "));
            Console.WriteLine();
        }

        foreach (var op in opsFaithful)
        {
            Console.WriteLine("############################# Faithful Group Action");
            var sdp = Group.SemiDirectProd(grN, op, grG);
            DisplayGroup.HeadSdp(sdp);
            Console.WriteLine(sdp.ElementsOrdersList().Glue(", "));
            Console.WriteLine();
        }

        var opFirst = opsGautN.FirstOrDefault(op => op?.Values.Distinct().Count() > 1, null);
        if (opsFaithful.Count == 0 && opsTransitives.Count == 0 && opFirst is not null)
        {
            Console.WriteLine("############################# First Group Action");
            var sdp = Group.SemiDirectProd(grN, opFirst, grG);
            DisplayGroup.HeadSdp(sdp);
            Console.WriteLine(sdp.ElementsOrdersList().Glue(", "));
            Console.WriteLine();
        }
    }

    static void CheckingPropertiesForSet<T>(ConcreteGroup<T> grG, int nbPoints) where T : struct, IElt<T>
    {
        var setX = new Symm(nbPoints);
        var opsGautN = Group.AllHomomorphisms(grG, setX);
        // TODO
    }

    private static Cn C2 => new Cn(2);
    private static Cn C3 => new Cn(3);
    private static Cn C4 => new Cn(4);
    private static Cn C5 => new Cn(5);
    private static Cn C6 => new Cn(6);
    private static Cn C7 => new Cn(7);
    private static Cn C8 => new Cn(8);
    private static ConcreteGroup<Ep2<ZnInt, ZnInt>> C3_2 => Product.Generate(C3, C3);
    private static ConcreteGroup<Ep2<ZnInt, ZnInt>> C4C2 => Product.Generate(C4, C2);

    public static void Examples()
    {
        // Neutral x = id is always fixed

        CheckPropertiesForSemiDirectProduct(C4, C2); // D8
        CheckPropertiesForSemiDirectProduct(C4, C2, rmPoints: C4[2]); // D8 transitive when x=0 and x=2 are removed
        CheckPropertiesForSemiDirectProduct(C5, C4); // Frobenius20
        CheckPropertiesForSemiDirectProduct(C7, C3, rmPoints: new[] { C7[1], C7[2], C7[4] });
        CheckPropertiesForSemiDirectProduct(C7, C6);
        CheckPropertiesForSemiDirectProduct(C3, C4); // No faithful action
        CheckPropertiesForSemiDirectProduct(C4, C4); // No faithful action
        CheckPropertiesForSemiDirectProduct(C4, C4, rmPoints: C4[2]); // No faithful action but transitive
        CheckPropertiesForSemiDirectProduct(C3, C8);
        CheckPropertiesForSemiDirectProduct(C4, C8, C4[2]);
    }
}