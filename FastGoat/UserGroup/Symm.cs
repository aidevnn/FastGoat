namespace FastGoat.UserGroup;

public class Symm : ConcreteGroup<Perm>
{
    public Symm(int n) : base($"Symm{n}", new Sn(n))
    {
        N = n;
    }

    public int N { get; }

    public IEnumerable<Perm[]> SubGroupsGenerators()
    {
        if (N >= 4)
        {
            // Direct Product
            yield return new[] { this[(1, 2)], this[(3, 4)] }; // V
            
            // Semi-direct Product
            yield return new[] { this[(1, 3)], this[(1, 2, 3, 4)] }; // D8
            
            // Symmetric 
            yield return new[] { this[(1, 2)], this[(1, 2, 3)] }; // S3
        }

        if (N >= 5)
        {
            // Direct Product
            yield return new[] { this[(1, 2)], this[(3, 4, 5)] }; // C6

            // Semi-direct Product
            yield return new[] { this[(2, 5), (3, 4)], this[(1, 2, 3, 4, 5)] }; // D10
            yield return new[] { this[(1, 2, 5)], this[(1, 2)], this[(3, 4)] }; // D12
            yield return new[] { this[(2, 3, 5, 4)], this[(1, 2, 3, 4, 5)] }; // C5 : C4

            // Symmetric or Alternate 
            yield return new[] { this[(1, 2), (3, 4)], this[(1, 2, 3)] }; // A4
            yield return new[] { this[(1, 2)], this[(1, 2, 3, 4)] }; // S4
        }

        if (N >= 6)
        {
            // Direct Product
            yield return new[] { this[(1, 2)], this[(3, 4, 5, 6)] }; // C8
            yield return new[] { this[(1, 2, 3)], this[(4, 5, 6)] }; // C3 x C3
            yield return new[] { this[(1, 2)], this[(3, 4), this[(5, 6)]] }; // C2 x C2 x C2

            // Semi-direct Product TODO 
        }

        if (N >= 7)
        {
            // Direct Product
            yield return new[] { this[(1, 2)], this[(3, 4, 5, 6, 7)] }; // C10
            yield return new[] { this[(1, 2, 3)], this[(4, 5, 6, 7)] }; // C12
            yield return new[] { this[(1, 2)], this[(3, 4), this[(5, 6, 7)]] }; // C2 x C2 x C3

            // Semi-direct Product TODO 
        }

        if (N >= 8)
        {
            // Direct Product
            yield return new[] { this[(1, 2)], this[(3, 4, 5, 6, 7, 8)] }; // C2 x C6
            yield return new[] { this[(1, 2, 3)], this[(4, 5, 6, 7, 8)] }; // C15
            yield return new[] { this[(1, 2, 3, 4)], this[(5, 6, 7, 8)] }; // C4 x C4
            yield return new[] { this[(1, 2)], this[(3, 4), this[(5, 6, 7, 8)]] }; // C2 x C2 x C4
            yield return new[] { this[(1, 2)], this[(3, 4), this[(5, 6)], this[(7, 8)]] }; // C2 x C2 x C2 x C2

            // Semi-direct Product TODO 
        }
    }
}