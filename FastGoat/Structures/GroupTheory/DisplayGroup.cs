using FastGoat.Structures.SetTheory;

namespace FastGoat.Structures.GroupTheory;


public class DisplayGroup<U> : DisplaySet<U> where U : struct, IElt<U>
{
    public DisplayGroup(SubGroup<U> group) : base(group)
    {
        Gr = group;
        SpecialChar(('@', Gr.Neutral));
    }

    SubGroup<U> Gr { get; }

    void ComputeOrders()
    {
        if (Gr.EltOrder)
            return;

        Gr.ComputeOrders();
    }

    public override void DisplayHead()
    {
        ComputeOrders();
        base.DisplayHead();
        Console.WriteLine("IsGroup      :{0,6}", Gr.IsGroup());
        Console.WriteLine("IsCommutative:{0,6}", Gr.IsCommutative());
        Console.WriteLine();
    }

    protected override string DisplayElement(U e, string name)
    {
        if (Gr.DisplayOrders)
            return string.Format("{0} = {1}", name, e);

        return string.Format("{0}{1,-4} = {2}", name, $"[{Gr.GetOrder(e)}]", e);
    }

    public override void DisplayElements()
    {
        ComputeOrders();
        base.DisplayElements();
    }

    public override void DisplayTable(char symb, Func<U, U, U> Op)
    {
        ComputeOrders();
        base.DisplayTable(symb, Op);
    }

}